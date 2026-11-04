//! Parameter estimation by weighted least squares, with profile likelihood.
//!
//! The objective is `chi2(p) = sum ((y_model(t) - y_obs) / sigma)^2`, minimized by Levenberg-Marquardt
//! in the logarithm of each positive parameter, with a Jacobian from forward differences over full
//! simulations. A profile fixes one parameter on a grid and re-optimizes the others; the points where
//! `chi2` exceeds its minimum by the 95 % quantile of chi-squared with one degree of freedom (3.84)
//! bound its confidence interval (Raue et al. 2009).

use super::compile::KineticModel;
use super::sim::{SimOptions, simulate};
use crate::error::{Error, Result};
use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};

/// One measurement.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Observation {
    /// Time of the measurement.
    pub time: f64,
    /// Symbol measured; `S#a` means the amount of species `S`.
    pub symbol: String,
    /// Measured value.
    pub value: f64,
    /// Standard deviation of the measurement.
    pub sigma: f64,
}

/// A parameter to estimate.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FitParameter {
    /// Symbol.
    pub name: String,
    /// Starting value, which must be positive.
    pub start: f64,
}

/// Result of a fit.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FitResult {
    /// Parameter names.
    pub parameters: Vec<String>,
    /// Estimates.
    pub estimates: Vec<f64>,
    /// Sum of squared weighted residuals at the optimum.
    pub chi2: f64,
    /// Number of observations.
    pub observations: usize,
    /// Akaike information criterion `chi2 + 2 k`.
    pub aic: f64,
    /// Levenberg-Marquardt iterations.
    pub iterations: usize,
    /// Standard errors of the logarithms from the Gauss-Newton Hessian, NaN if singular.
    pub log_standard_errors: Vec<f64>,
}

struct Problem<'a> {
    base: &'a KineticModel,
    names: Vec<String>,
    data: &'a [Observation],
    times: Vec<f64>,
    fixed: Vec<(String, f64)>,
    options: SimOptions,
}

impl Problem<'_> {
    fn residuals(&self, log_p: &[f64]) -> Option<Vec<f64>> {
        let mut m = self.base.clone();
        for (name, v) in &self.fixed {
            m.set_value(name, *v).ok()?;
        }
        for (name, lp) in self.names.iter().zip(log_p) {
            m.set_value(name, lp.exp()).ok()?;
        }
        let course = simulate(&m, &self.times, &self.options).ok()?;
        let mut out = Vec::with_capacity(self.data.len());
        for obs in self.data {
            let k = self.times.iter().position(|t| *t == obs.time)?;
            let slot = m.slot(&obs.symbol).ok()?;
            let v = course.values[k][slot];
            if !v.is_finite() {
                return None;
            }
            out.push((v - obs.value) / obs.sigma);
        }
        Some(out)
    }

    fn chi2(&self, log_p: &[f64]) -> Option<f64> {
        self.residuals(log_p).map(|r| r.iter().map(|v| v * v).sum())
    }

    fn solve(&self, start: Vec<f64>, max_iterations: usize) -> Result<(Vec<f64>, f64, usize, DMatrix<f64>)> {
        let mut p = start;
        let mut r = self
            .residuals(&p)
            .ok_or_else(|| Error::Kinetic("the model cannot be simulated at the starting values".into()))?;
        let mut cost: f64 = r.iter().map(|v| v * v).sum();
        let mut lambda = 1e-3;
        let k = p.len();
        let mut jtj = DMatrix::zeros(k, k);
        let mut iterations = 0;
        for _ in 0..max_iterations {
            iterations += 1;
            let mut jac = DMatrix::zeros(r.len(), k);
            for c in 0..k {
                let mut q = p.clone();
                let h = 1e-5 * (1.0 + p[c].abs());
                q[c] += h;
                let rq = self.residuals(&q).ok_or_else(|| Error::Kinetic("simulation failed during the fit".into()))?;
                for i in 0..r.len() {
                    jac[(i, c)] = (rq[i] - r[i]) / h;
                }
            }
            jtj = jac.transpose() * &jac;
            let gradient = jac.transpose() * DVector::from_vec(r.clone());
            let mut improved = false;
            for _ in 0..30 {
                let mut a = jtj.clone();
                for d in 0..k {
                    a[(d, d)] += lambda * jtj[(d, d)].max(1e-12);
                }
                if let Some(step) = a.lu().solve(&(-&gradient)) {
                    let q: Vec<f64> = p.iter().zip(step.iter()).map(|(x, s)| x + s).collect();
                    if let Some(rq) = self.residuals(&q) {
                        let new_cost: f64 = rq.iter().map(|v| v * v).sum();
                        if new_cost < cost {
                            let drop = cost - new_cost;
                            p = q;
                            r = rq;
                            cost = new_cost;
                            lambda = (lambda * 0.3).max(1e-12);
                            improved = true;
                            if drop < 1e-10 * (1.0 + cost) {
                                return Ok((p, cost, iterations, jtj));
                            }
                            break;
                        }
                    }
                }
                lambda *= 4.0;
                if lambda > 1e12 {
                    break;
                }
            }
            if !improved {
                break;
            }
        }
        Ok((p, cost, iterations, jtj))
    }
}

fn observation_times(data: &[Observation]) -> Vec<f64> {
    let mut times: Vec<f64> = std::iter::once(0.0).chain(data.iter().map(|o| o.time)).collect();
    times.sort_by(f64::total_cmp);
    times.dedup();
    times
}

/// Fits `parameters` to `data`.
pub fn fit(
    model: &KineticModel,
    parameters: &[FitParameter],
    data: &[Observation],
    options: &SimOptions,
) -> Result<FitResult> {
    if parameters.iter().any(|p| p.start <= 0.0) {
        return Err(Error::Invalid("starting values must be positive".into()));
    }
    let problem = Problem {
        base: model,
        names: parameters.iter().map(|p| p.name.clone()).collect(),
        data,
        times: observation_times(data),
        fixed: Vec::new(),
        options: *options,
    };
    let start: Vec<f64> = parameters.iter().map(|p| p.start.ln()).collect();
    let (p, chi2, iterations, jtj) = problem.solve(start, 200)?;
    let errors = match jtj.lu().try_inverse() {
        Some(inv) => (0..p.len()).map(|i| inv[(i, i)].abs().sqrt()).collect(),
        None => vec![f64::NAN; p.len()],
    };
    Ok(FitResult {
        parameters: problem.names.clone(),
        estimates: p.iter().map(|v| v.exp()).collect(),
        chi2,
        observations: data.len(),
        aic: chi2 + 2.0 * parameters.len() as f64,
        iterations,
        log_standard_errors: errors,
    })
}

/// A likelihood profile of one parameter.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Profile {
    /// Profiled parameter.
    pub parameter: String,
    /// Grid values.
    pub values: Vec<f64>,
    /// `chi2` minimized over the other parameters at each grid value.
    pub chi2: Vec<f64>,
    /// Lower limit of the 95 % interval, when the profile crosses the threshold on that side.
    pub lower: Option<f64>,
    /// Upper limit of the 95 % interval, when the profile crosses the threshold on that side.
    pub upper: Option<f64>,
}

const CHI2_95_1DOF: f64 = 3.841_458_820_694_124;

/// Profiles `name` over `factor` either side of `estimate`, on `points` geometric grid values.
#[allow(clippy::too_many_arguments)]
pub fn profile(
    model: &KineticModel,
    parameters: &[FitParameter],
    data: &[Observation],
    fit_result: &FitResult,
    name: &str,
    factor: f64,
    points: usize,
    options: &SimOptions,
) -> Result<Profile> {
    let at = fit_result
        .parameters
        .iter()
        .position(|p| p == name)
        .ok_or_else(|| Error::Unknown { kind: "parameter", id: name.into() })?;
    let estimate = fit_result.estimates[at];
    let others: Vec<FitParameter> = parameters
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != at)
        .map(|(i, p)| FitParameter { name: p.name.clone(), start: fit_result.estimates[i] })
        .collect();
    let times = observation_times(data);
    let grid: Vec<f64> = (0..points)
        .map(|k| estimate * factor.powf(2.0 * k as f64 / (points - 1).max(1) as f64 - 1.0))
        .collect();
    let chi2: Vec<f64> = crate::parallel::map_chunks(&grid, 1, |slice| {
        let problem = Problem {
            base: model,
            names: others.iter().map(|p| p.name.clone()).collect(),
            data,
            times: times.clone(),
            fixed: vec![(name.to_string(), slice[0])],
            options: *options,
        };
        let start: Vec<f64> = others.iter().map(|p| p.start.ln()).collect();
        if start.is_empty() {
            problem.chi2(&[]).unwrap_or(f64::INFINITY)
        } else {
            problem.solve(start, 50).map_or(f64::INFINITY, |r| r.1)
        }
    });
    let threshold = fit_result.chi2 + CHI2_95_1DOF;
    let centre = grid.iter().position(|v| *v >= estimate).unwrap_or(points / 2);
    let crossing = |range: Vec<usize>| -> Option<f64> {
        let mut prev: Option<usize> = None;
        for k in range {
            if chi2[k] > threshold {
                return Some(match prev {
                    Some(p) => {
                        let w = (threshold - chi2[p]) / (chi2[k] - chi2[p]);
                        (grid[p].ln() + w * (grid[k].ln() - grid[p].ln())).exp()
                    }
                    None => grid[k],
                });
            }
            prev = Some(k);
        }
        None
    };
    let lower = crossing((0..=centre.min(points - 1)).rev().collect());
    let upper = crossing((centre.min(points - 1)..points).collect());
    Ok(Profile { parameter: name.to_string(), values: grid, chi2, lower, upper })
}
