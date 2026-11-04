//! Metabolic control analysis at a steady state.
//!
//! With `N = L Nr` and unscaled elasticities `E = dv/dy` from symbolic differentiation of the rate
//! laws, the concentration and flux control matrices are (Reder 1988):
//! `Cs = -L (Nr E L)^-1 Nr` and `Cj = I + E Cs`. Scaled coefficients multiply by `v_j / y_i` and
//! `v_j / v_i`. The summation theorems read `sum_j Cj_ij = 1` and `sum_j Cs_ij = 0`, the connectivity
//! theorems `Cj E L = 0` and `Cs E L = -L`. [`ControlAnalysis::theorems`] reports their residuals in
//! the scaled form in which the theorems are stated, and the unscaled residuals beside them.

use super::compile::{KineticModel, NoDelay, compile_expr, eval};
use super::diff::{Var, diff};
use super::expr::Expr;
use super::steady::{Moieties, SteadyState, moieties};
use crate::error::{Error, Result};
use crate::linalg::nullspace;
use nalgebra::DMatrix;
use serde::Serialize;

/// Residuals of the control theorems.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct Theorems {
    /// Largest deviation in `sum_j Cj_ij = 1`, over reactions with non-zero flux, scaled.
    pub flux_summation: f64,
    /// Largest deviation in `sum_j Cs_ij = 0` over species with non-zero amount, scaled.
    pub concentration_summation: f64,
    /// Largest scaled entry of `Cj E L`: the connectivity theorem `sum_j C^J_j e^j_k = 0`.
    pub flux_connectivity: f64,
    /// Largest scaled entry of `Cs E L + L`: the connectivity theorem `sum_j C^S_j e^j_k = -delta`.
    pub concentration_connectivity: f64,
    /// Largest entry of `Cj E L` as computed, which carries the units of the elasticities.
    pub unscaled_flux_connectivity: f64,
    /// Largest entry of `Cs E L + L` as computed.
    pub unscaled_concentration_connectivity: f64,
    /// Largest entry of `Cj K - K` for the steady-state flux modes `K`.
    pub flux_kernel: f64,
    /// Largest entry of `Cs K`.
    pub concentration_kernel: f64,
}

/// Control analysis at a steady state.
#[derive(Clone, Debug)]
pub struct ControlAnalysis {
    /// Reaction identifiers.
    pub reactions: Vec<String>,
    /// State names.
    pub states: Vec<String>,
    /// Fluxes.
    pub fluxes: Vec<f64>,
    /// State values.
    pub y: Vec<f64>,
    /// Unscaled elasticities, `reactions x states`.
    pub elasticities: DMatrix<f64>,
    /// Scaled elasticities.
    pub scaled_elasticities: DMatrix<f64>,
    /// Unscaled flux control, `reactions x reactions`.
    pub flux_control: DMatrix<f64>,
    /// Scaled flux control.
    pub scaled_flux_control: DMatrix<f64>,
    /// Unscaled concentration control, `states x reactions`.
    pub concentration_control: DMatrix<f64>,
    /// Scaled concentration control.
    pub scaled_concentration_control: DMatrix<f64>,
    /// Residuals of the theorems.
    pub theorems: Theorems,
    /// The factorization used.
    pub moieties: Moieties,
}

fn max_abs(m: &DMatrix<f64>) -> f64 {
    m.iter().fold(0.0f64, |a, v| a.max(v.abs()))
}

/// Solves `A X = B` with row and column equilibration and two steps of iterative refinement. Enzymes
/// near equilibrium make `A` badly scaled, by many orders of magnitude between rows.
fn solve_refined(a: &DMatrix<f64>, b: &DMatrix<f64>) -> Option<DMatrix<f64>> {
    let n = a.nrows();
    let rows: Vec<f64> =
        (0..n).map(|i| 1.0 / a.row(i).iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-300)).collect();
    let cols: Vec<f64> =
        (0..n).map(|j| 1.0 / (0..n).fold(0.0f64, |m, i| m.max((a[(i, j)] * rows[i]).abs())).max(1e-300)).collect();
    let scaled = DMatrix::from_fn(n, n, |i, j| rows[i] * a[(i, j)] * cols[j]);
    let lu = scaled.lu();
    let solve_one = |rhs: &DMatrix<f64>| -> Option<DMatrix<f64>> {
        let scaled_rhs = DMatrix::from_fn(n, rhs.ncols(), |i, j| rows[i] * rhs[(i, j)]);
        let y = lu.solve(&scaled_rhs)?;
        Some(DMatrix::from_fn(n, rhs.ncols(), |i, j| cols[i] * y[(i, j)]))
    };
    let mut x = solve_one(b)?;
    for _ in 0..2 {
        let residual = b - a * &x;
        x += solve_one(&residual)?;
    }
    Some(x)
}

/// Unscaled elasticity matrix at given slot values, by symbolic differentiation.
pub fn elasticities(model: &KineticModel, values: &[f64]) -> Result<DMatrix<f64>> {
    if model.flat_rates.len() != model.reactions.len() {
        return Err(Error::Kinetic("rate laws were too large to differentiate".into()));
    }
    let mut e = DMatrix::zeros(model.reactions.len(), model.n_states());
    for (j, law) in model.flat_rates.iter().enumerate() {
        let mut syms = Vec::new();
        law.symbols(&mut syms);
        syms.sort();
        syms.dedup();
        for name in syms {
            if let Some(k) = model.state_names.iter().position(|n| *n == name) {
                let d = diff(law, Var::Sym(&name));
                if !matches!(d, Expr::Num(v) if v == 0.0) {
                    let c = compile_expr(&d, &model.index)?;
                    e[(j, k)] = eval(&c, values, 0.0, &NoDelay);
                }
            }
        }
    }
    Ok(e)
}

/// Partial derivatives of each rate with respect to constant symbols, `reactions x parameters`.
pub fn parameter_elasticities(model: &KineticModel, values: &[f64], parameters: &[String]) -> Result<DMatrix<f64>> {
    let mut e = DMatrix::zeros(model.reactions.len(), parameters.len());
    for (j, law) in model.flat_rates.iter().enumerate() {
        for (k, name) in parameters.iter().enumerate() {
            if law.mentions(name) {
                let d = diff(law, Var::Sym(name));
                if !matches!(d, Expr::Num(v) if v == 0.0) {
                    let c = compile_expr(&d, &model.index)?;
                    e[(j, k)] = eval(&c, values, 0.0, &NoDelay);
                }
            }
        }
    }
    Ok(e)
}

/// Runs control analysis at `state`.
pub fn control_analysis(model: &KineticModel, state: &SteadyState) -> Result<ControlAnalysis> {
    let m = moieties(model)?;
    let e = elasticities(model, &state.values)?;
    let reduced = &m.nr * &e * &m.link;
    let solved = solve_refined(&reduced, &m.nr)
        .ok_or_else(|| Error::Kinetic("the reduced Jacobian is singular at this steady state".into()))?;
    let cs = -(&m.link * solved);
    let n_r = model.reactions.len();
    let n_s = model.n_states();
    let cj = DMatrix::identity(n_r, n_r) + &e * &cs;

    let v = &state.fluxes;
    let y = &state.y;
    let mut scaled_e = DMatrix::zeros(n_r, n_s);
    for j in 0..n_r {
        for k in 0..n_s {
            scaled_e[(j, k)] = if v[j] != 0.0 { e[(j, k)] * y[k] / v[j] } else { f64::NAN };
        }
    }
    let mut scaled_cj = DMatrix::zeros(n_r, n_r);
    for i in 0..n_r {
        for j in 0..n_r {
            scaled_cj[(i, j)] = if v[i] != 0.0 { cj[(i, j)] * v[j] / v[i] } else { f64::NAN };
        }
    }
    let mut scaled_cs = DMatrix::zeros(n_s, n_r);
    for i in 0..n_s {
        for j in 0..n_r {
            scaled_cs[(i, j)] = if y[i] != 0.0 { cs[(i, j)] * v[j] / y[i] } else { f64::NAN };
        }
    }

    // A flux that is zero to rounding error (a reaction at equilibrium) has no scaled coefficients.
    let flux_scale = v.iter().fold(0.0f64, |a, x| a.max(x.abs()));
    let carries_flux = |i: usize| v[i].abs() > 1e-9 * flux_scale;

    // Theorem residuals, in the scaled form in which the theorems are stated.
    let mut flux_summation = 0.0f64;
    for i in 0..n_r {
        if carries_flux(i) {
            let sum: f64 = (0..n_r).map(|j| scaled_cj[(i, j)]).sum();
            flux_summation = flux_summation.max((sum - 1.0).abs());
        }
    }
    let mut concentration_summation = 0.0f64;
    let dynamic: Vec<usize> = (0..n_s).filter(|&i| (0..n_r).any(|j| m.n[(i, j)] != 0.0)).collect();
    for &i in &dynamic {
        if y[i] != 0.0 {
            let sum: f64 = (0..n_r).map(|j| scaled_cs[(i, j)]).sum();
            concentration_summation = concentration_summation.max(sum.abs());
        }
    }
    let kernel = nullspace(&m.n, 1e-10);
    let flux_conn = &cj * &e * &m.link;
    let conc_conn = &cs * &e * &m.link + &m.link;
    // Scaled connectivity: `(y_k / v_i) (Cj E L)_ik` and `(y_k / y_i) (Cs E L + L)_ik` over the
    // independent species `k`.
    let mut flux_connectivity = 0.0f64;
    for i in 0..n_r {
        if carries_flux(i) {
            for (col, &k) in m.independent.iter().enumerate() {
                flux_connectivity = flux_connectivity.max((y[k] / v[i] * flux_conn[(i, col)]).abs());
            }
        }
    }
    let mut concentration_connectivity = 0.0f64;
    for &i in &dynamic {
        if y[i] != 0.0 {
            for (col, &k) in m.independent.iter().enumerate() {
                concentration_connectivity = concentration_connectivity.max((y[k] / y[i] * conc_conn[(i, col)]).abs());
            }
        }
    }
    let theorems = Theorems {
        flux_summation,
        concentration_summation,
        flux_connectivity,
        concentration_connectivity,
        unscaled_flux_connectivity: max_abs(&flux_conn),
        unscaled_concentration_connectivity: max_abs(&conc_conn),
        flux_kernel: max_abs(&(&cj * &kernel - &kernel)),
        concentration_kernel: max_abs(&(&cs * &kernel)),
    };
    Ok(ControlAnalysis {
        reactions: model.reactions.iter().map(|r| r.id.clone()).collect(),
        states: model.state_names.clone(),
        fluxes: v.clone(),
        y: y.clone(),
        elasticities: e,
        scaled_elasticities: scaled_e,
        flux_control: cj,
        scaled_flux_control: scaled_cj,
        concentration_control: cs,
        scaled_concentration_control: scaled_cs,
        theorems,
        moieties: m,
    })
}
