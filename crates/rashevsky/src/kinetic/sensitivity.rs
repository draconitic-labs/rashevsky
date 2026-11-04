//! Parameter sensitivities: local steady-state coefficients and variance-based global indices.

use super::compile::KineticModel;
use super::mca::{ControlAnalysis, parameter_elasticities};
use super::steady::SteadyState;
use crate::error::{Error, Result};
use nalgebra::DMatrix;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};

/// Local steady-state sensitivities: `dJ/dp = Cj dv/dp` and `dy/dp = Cs dv/dp`.
#[derive(Clone, Debug)]
pub struct LocalSensitivity {
    /// Parameter names.
    pub parameters: Vec<String>,
    /// Flux sensitivities, `reactions x parameters`.
    pub flux: DMatrix<f64>,
    /// State sensitivities, `states x parameters`.
    pub state: DMatrix<f64>,
}

/// Local steady-state sensitivities with respect to constant symbols.
pub fn local_steady(
    model: &KineticModel,
    state: &SteadyState,
    control: &ControlAnalysis,
    parameters: &[String],
) -> Result<LocalSensitivity> {
    let dv = parameter_elasticities(model, &state.values, parameters)?;
    Ok(LocalSensitivity {
        parameters: parameters.to_vec(),
        flux: &control.flux_control * &dv,
        state: &control.concentration_control * &dv,
    })
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
    /// Samples that failed to evaluate and were replaced by resampling.
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
    if n < 2 {
        return Err(Error::Kinetic("too few Sobol samples evaluated".into()));
    }
    let mean: f64 = rows.iter().map(|r| 0.5 * (r.0 + r.1)).sum::<f64>() / n as f64;
    let variance: f64 =
        rows.iter().map(|r| 0.5 * ((r.0 - mean).powi(2) + (r.1 - mean).powi(2))).sum::<f64>() / n as f64;
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
