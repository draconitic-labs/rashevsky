//! Parameter sensitivities: local steady-state coefficients and variance-based global indices.

use super::compile::{KineticModel, NoDelay};
use super::mca::ControlAnalysis;
use super::sim::{SimOptions, simulate};
use super::steady::{SteadyOptions, SteadyState, steady_state};
use crate::error::{Error, Result};
use nalgebra::{DMatrix, DVector};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};

/// Local steady-state sensitivities. A parameter moves the rates of change at the steady state by
/// `b = N dv/dp + (dN/dp) v`, and it can move the start: with `d0 = dy0/dp - L dy0_ind/dp`, the part of that
/// move that the dynamics cannot undo, such as a conserved total, the independent amounts respond by
/// `ds/dp = -(Nr E L)^-1 (b + Nr E d0)`, which gives `dy/dp = d0 + L ds/dp` and `dJ/dp = dv/dp + E dy/dp`.
/// The partial derivatives at the frozen state come from a
/// fourth-order central difference through [`KineticModel::set_value`], so a parameter counts everywhere
/// it reaches: laws, stoichiometries, compartment sizes and other values assigned from it.
#[derive(Clone, Debug)]
pub struct LocalSensitivity {
    /// Parameter names.
    pub parameters: Vec<String>,
    /// Flux sensitivities, `reactions x parameters`.
    pub flux: DMatrix<f64>,
    /// State sensitivities, `states x parameters`.
    pub state: DMatrix<f64>,
}

/// Rates of change and reaction rates at the frozen `state` with one parameter overridden.
fn frozen_rates(
    model: &KineticModel,
    state: &SteadyState,
    name: &str,
    value: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let mut m = model.clone();
    m.set_value(name, value)?;
    let mut vals = m.initial.clone();
    m.eval_state(state.time, &state.y, &mut vals, &NoDelay);
    let mut dy = vec![0.0; state.y.len()];
    m.rhs_from_vals(state.time, &vals, &mut dy, &NoDelay);
    let rates = m.reactions.iter().map(|r| vals[r.rate]).collect();
    Ok((dy, rates, m.initial_state()))
}

/// Local steady-state sensitivities with respect to constant symbols.
pub fn local_steady(
    model: &KineticModel,
    state: &SteadyState,
    control: &ControlAnalysis,
    parameters: &[String],
) -> Result<LocalSensitivity> {
    let varying = model.varying_slots();
    let moieties = &control.moieties;
    let ind = &moieties.independent;
    let reduced = &moieties.nr * &control.elasticities * &moieties.link;
    let lu = reduced.lu();
    let (n_states, n_reactions, k) = (state.y.len(), model.reactions.len(), parameters.len());
    let mut flux = DMatrix::zeros(n_reactions, k);
    let mut response = DMatrix::zeros(n_states, k);
    for (column, name) in parameters.iter().enumerate() {
        let slot = model.slot(name)?;
        if varying.contains(&slot) {
            return Err(Error::Invalid(format!("`{name}` varies during a simulation, so it has no local sensitivity")));
        }
        let p = model.initial[slot];
        let h = 1e-3 * if p == 0.0 { 1.0 } else { p.abs() };
        let at = |offset: f64| frozen_rates(model, state, name, p + offset * h);
        let (a, b, c, d) = (at(-2.0)?, at(-1.0)?, at(1.0)?, at(2.0)?);
        let stencil = |w: f64, x: f64, y: f64, z: f64| (w - 8.0 * x + 8.0 * y - z) / (12.0 * h);
        let dv = DVector::from_iterator(n_reactions, (0..n_reactions).map(|j| stencil(a.1[j], b.1[j], c.1[j], d.1[j])));
        let dy0 = DVector::from_iterator(n_states, (0..n_states).map(|i| stencil(a.2[i], b.2[i], c.2[i], d.2[i])));
        let dy0_ind = DVector::from_iterator(ind.len(), ind.iter().map(|&i| dy0[i]));
        let unreachable = &dy0 - &moieties.link * dy0_ind;
        let through_rates = &moieties.nr * (&control.elasticities * &unreachable);
        let rhs = DVector::from_iterator(
            ind.len(),
            ind.iter().enumerate().map(|(row, &i)| -stencil(a.0[i], b.0[i], c.0[i], d.0[i]) - through_rates[row]),
        );
        let ds = if ind.is_empty() {
            DVector::zeros(0)
        } else {
            lu.solve(&rhs)
                .ok_or_else(|| Error::Kinetic("the reduced Jacobian is singular at this steady state".into()))?
        };
        let dy = &moieties.link * ds + unreachable;
        let dj = dv + &control.elasticities * &dy;
        flux.column_mut(column).copy_from(&dj);
        response.column_mut(column).copy_from(&dy);
    }
    Ok(LocalSensitivity { parameters: parameters.to_vec(), flux, state: response })
}

/// A parameter's range for global analysis.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Range {
    /// Name.
    pub name: String,
    /// Lower limit.
    pub low: f64,
    /// Upper limit.
    pub high: f64,
    /// Sample uniformly in the logarithm rather than the value.
    #[serde(default)]
    pub log: bool,
}

/// First-order and total Sobol indices of a scalar output.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SobolIndices {
    /// Parameter names.
    pub parameters: Vec<String>,
    /// First-order index per parameter.
    pub first: Vec<f64>,
    /// Total-order index per parameter.
    pub total: Vec<f64>,
    /// Output evaluations used.
    pub evaluations: usize,
    /// Samples that failed to evaluate and were left out. A failure region biases the indices, so more
    /// than 5 % of samples failing is an error.
    pub failures: usize,
}

fn draw(ranges: &[Range], u: &[f64]) -> Vec<f64> {
    ranges
        .iter()
        .zip(u)
        .map(
            |(r, &v)| {
                if r.log { (r.low.ln() + v * (r.high.ln() - r.low.ln())).exp() } else { r.low + v * (r.high - r.low) }
            },
        )
        .collect()
}

/// Sobol indices by the estimators of Saltelli et al. (2010): `S_i` from their equation (b) and
/// `T_i` from Jansen's formula, with `n (d + 2)` evaluations of `f`.
pub fn sobol(
    ranges: &[Range],
    samples: usize,
    seed: u64,
    f: &(dyn Fn(&[f64]) -> Option<f64> + Sync),
) -> Result<SobolIndices> {
    let d = ranges.len();
    if d == 0 || samples < 2 {
        return Err(Error::Invalid("Sobol analysis needs parameters and at least two samples".into()));
    }
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut a = vec![vec![0.0; d]; samples];
    let mut b = vec![vec![0.0; d]; samples];
    for row in a.iter_mut().chain(b.iter_mut()) {
        for v in row.iter_mut() {
            *v = rng.random();
        }
    }
    // Evaluate in parallel; a sample with a failed evaluation in any arm is dropped from every arm.
    let idx: Vec<usize> = (0..samples).collect();
    type Row = Option<(f64, f64, Vec<f64>)>;
    let parts: Vec<Vec<Row>> = crate::parallel::map_chunks(&idx, crate::parallel::chunk_size(samples), |slice| {
        slice
            .iter()
            .map(|&k| {
                let fa = f(&draw(ranges, &a[k]))?;
                let fb = f(&draw(ranges, &b[k]))?;
                let mut fab = Vec::with_capacity(d);
                for i in 0..d {
                    let mut mix = a[k].clone();
                    mix[i] = b[k][i];
                    fab.push(f(&draw(ranges, &mix))?);
                }
                Some((fa, fb, fab))
            })
            .collect()
    });
    let rows: Vec<(f64, f64, Vec<f64>)> = parts.into_iter().flatten().flatten().collect();
    let n = rows.len();
    if n < 2 || (samples - n) * 20 > samples {
        return Err(Error::Kinetic(format!("{} of {samples} Sobol samples could not be evaluated", samples - n)));
    }
    let mean: f64 = rows.iter().map(|r| 0.5 * (r.0 + r.1)).sum::<f64>() / n as f64;
    let variance: f64 =
        rows.iter().map(|r| 0.5 * ((r.0 - mean).powi(2) + (r.1 - mean).powi(2))).sum::<f64>() / n as f64;
    // An output that does not respond has variance of the size of rounding, and ratios of that mean nothing.
    let floor = (1e-10 * mean.abs()).powi(2);
    if variance <= floor {
        return Ok(SobolIndices {
            parameters: ranges.iter().map(|r| r.name.clone()).collect(),
            first: vec![0.0; d],
            total: vec![0.0; d],
            evaluations: n * (d + 2),
            failures: samples - n,
        });
    }
    let mut first = vec![0.0; d];
    let mut total = vec![0.0; d];
    for i in 0..d {
        let mut s = 0.0;
        let mut t = 0.0;
        for (fa, fb, fab) in &rows {
            s += fb * (fab[i] - fa);
            t += (fa - fab[i]).powi(2);
        }
        first[i] = s / n as f64 / variance;
        total[i] = 0.5 * t / n as f64 / variance;
    }
    Ok(SobolIndices {
        parameters: ranges.iter().map(|r| r.name.clone()).collect(),
        first,
        total,
        evaluations: n * (d + 2),
        failures: samples - n,
    })
}

/// The scalar a global analysis of a kinetic model reads from each parameter sample.
#[derive(Clone, Debug)]
pub enum Output {
    /// A symbol at the steady state reached after integrating for `integrate`.
    Steady {
        /// Species, reaction or parameter.
        symbol: String,
        /// Integration before Newton iteration.
        integrate: f64,
    },
    /// A symbol at a time.
    AtTime {
        /// Species, reaction or parameter.
        symbol: String,
        /// Time of the reading.
        time: f64,
    },
}

/// Sobol indices of one output of a kinetic model over parameter ranges.
pub fn sobol_model(
    model: &KineticModel,
    ranges: &[Range],
    samples: usize,
    seed: u64,
    output: &Output,
) -> Result<SobolIndices> {
    for r in ranges {
        model.slot(&r.name)?;
        if r.low >= r.high || (r.log && r.low <= 0.0) {
            return Err(Error::Invalid(format!("range of `{}` is empty or not positive on a log scale", r.name)));
        }
    }
    let symbol = match output {
        Output::Steady { symbol, .. } | Output::AtTime { symbol, .. } => symbol,
    };
    let slot = model.slot(symbol)?;
    let evaluate = |values: &[f64]| -> Option<f64> {
        let mut m = model.clone();
        for (r, &v) in ranges.iter().zip(values) {
            m.set_value(&r.name, v).ok()?;
        }
        let reading = match output {
            Output::Steady { integrate, .. } => {
                let options = SteadyOptions { integrate_first: *integrate, ..SteadyOptions::default() };
                // A root that is not an attractor is not what this parameter set does.
                let state = steady_state(&m, &options).ok()?;
                if !state.stable {
                    return None;
                }
                state.values[slot]
            }
            Output::AtTime { time, .. } => {
                let run = simulate(&m, &[0.0, *time], &SimOptions::default()).ok()?;
                run.values.last()?[slot]
            }
        };
        reading.is_finite().then_some(reading)
    };
    sobol(ranges, samples, seed, &evaluate)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The Ishigami function has closed-form Sobol indices (Ishigami and Homma 1990).
    #[test]
    fn ishigami_indices_match_their_closed_form() {
        let pi = std::f64::consts::PI;
        let ranges: Vec<Range> =
            (0..3).map(|i| Range { name: format!("x{i}"), low: -pi, high: pi, log: false }).collect();
        let (a, b) = (7.0f64, 0.1f64);
        let result = sobol(&ranges, 60_000, 42, &|x: &[f64]| {
            Some(x[0].sin() + a * x[1].sin().powi(2) + b * x[2].powi(4) * x[0].sin())
        })
        .unwrap();
        let variance = 0.5 + a * a / 8.0 + b * pi.powi(4) / 5.0 + b * b * pi.powi(8) / 18.0;
        let v1 = 0.5 * (1.0 + b * pi.powi(4) / 5.0).powi(2);
        let v2 = a * a / 8.0;
        let v13 = b * b * pi.powi(8) * (1.0 / 18.0 - 1.0 / 50.0);
        let expected_first = [v1 / variance, v2 / variance, 0.0];
        let expected_total = [(v1 + v13) / variance, v2 / variance, v13 / variance];
        for i in 0..3 {
            assert!((result.first[i] - expected_first[i]).abs() < 0.02, "S{i}: {:?}", result.first);
            assert!((result.total[i] - expected_total[i]).abs() < 0.02, "T{i}: {:?}", result.total);
        }
    }
}
