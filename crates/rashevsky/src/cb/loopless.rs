//! Thermodynamically infeasible cycles.
//!
//! [`loopless_fba`] is the mixed-integer ll-FBA of Schellenberger, Lewis and Palsson (2011): it adds
//! a Gibbs-energy-like potential per internal reaction whose sign must oppose the flux, with
//! `N_int' G = 0` over the null space of the internal stoichiometry. [`cycle_free`] is the
//! CycleFreeFlux linear program of Desouki et al. (2015): it removes loops from a given flux vector
//! while holding exchanges fixed and keeping each internal flux between zero and its original value.

use super::fba::FluxSolution;
use crate::error::{Error, Result};
use crate::linalg::nullspace;
use crate::lp::{self, Backend, Problem, Sense, Status};
use crate::model::Model;
use nalgebra::DMatrix;

/// Reactions that are not boundary reactions, by index.
pub fn internal_reactions(model: &Model) -> Vec<usize> {
    model.reactions.iter().enumerate().filter(|(_, r)| !r.is_boundary()).map(|(j, _)| j).collect()
}

/// Adds the ll-FBA constraints to `p`. `big_m` bounds the fluxes and potentials used in the big-M rows.
fn add_loopless_constraints(model: &Model, p: &mut Problem, big_m: f64) {
    let internal = internal_reactions(model);
    let mut s = DMatrix::<f64>::zeros(model.metabolites.len(), internal.len());
    for (k, &j) in internal.iter().enumerate() {
        for &(i, c) in &model.reactions[j].stoichiometry {
            if !model.metabolites[i].boundary {
                s[(i, k)] = c;
            }
        }
    }
    let null = nullspace(&s, 1e-10);
    let mut g_cols = Vec::with_capacity(internal.len());
    for &j in &internal {
        let a = p.add_col(0.0, 0.0, 1.0, true);
        let g = p.add_col(0.0, -big_m, big_m, false);
        g_cols.push(g);
        // a = 1 -> v >= 0 and G <= -1; a = 0 -> v <= 0 and G >= 1.
        let (lb, ub) = (model.reactions[j].lb, model.reactions[j].ub);
        p.add_row(f64::NEG_INFINITY, 0.0, vec![(j, 1.0), (a, -ub.max(0.0))]);
        p.add_row(lb.min(0.0), f64::INFINITY, vec![(j, 1.0), (a, lb.min(0.0))]);
        p.add_row(f64::NEG_INFINITY, big_m - 1.0, vec![(g, 1.0), (a, big_m)]);
        p.add_row(1.0, f64::INFINITY, vec![(g, 1.0), (a, big_m + 1.0)]);
    }
    for k in 0..null.ncols() {
        let coefs: Vec<(usize, f64)> =
            (0..internal.len()).filter(|&i| null[(i, k)].abs() > 1e-12).map(|i| (g_cols[i], null[(i, k)])).collect();
        if !coefs.is_empty() {
            p.add_row(0.0, 0.0, coefs);
        }
    }
}

/// ll-FBA with the default backend. `big_m` is 1000 in Schellenberger et al.
pub fn loopless_fba(model: &Model, big_m: f64) -> Result<FluxSolution> {
    loopless_fba_with(model, big_m, lp::default_backend().as_ref())
}

/// ll-FBA with a chosen backend, which must solve mixed-integer programs.
pub fn loopless_fba_with(model: &Model, big_m: f64, backend: &dyn Backend) -> Result<FluxSolution> {
    let mut p = model.to_lp();
    add_loopless_constraints(model, &mut p, big_m);
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    Ok(FluxSolution {
        status: s.status,
        objective: s.objective,
        fluxes: s.x[..model.reactions.len()].to_vec(),
        shadow_prices: Vec::new(),
        reduced_costs: Vec::new(),
    })
}

/// CycleFreeFlux with the default backend.
pub fn cycle_free(model: &Model, fluxes: &[f64]) -> Result<Vec<f64>> {
    cycle_free_with(model, fluxes, lp::default_backend().as_ref())
}

/// CycleFreeFlux with a chosen backend.
pub fn cycle_free_with(model: &Model, fluxes: &[f64], backend: &dyn Backend) -> Result<Vec<f64>> {
    if fluxes.len() != model.reactions.len() {
        return Err(Error::Invalid("flux vector length differs from reaction count".into()));
    }
    let mut p = model.to_lp();
    p.sense = Some(Sense::Minimize);
    for (j, r) in model.reactions.iter().enumerate() {
        let v = fluxes[j];
        if r.is_boundary() {
            p.col_lb[j] = v;
            p.col_ub[j] = v;
            p.obj[j] = 0.0;
        } else if v >= 0.0 {
            p.col_lb[j] = 0.0;
            p.col_ub[j] = v;
            p.obj[j] = 1.0;
        } else {
            p.col_lb[j] = v;
            p.col_ub[j] = 0.0;
            p.obj[j] = -1.0;
        }
    }
    // Reactions in the objective keep their value, so the objective is unchanged.
    for &(j, _) in &model.objective {
        p.col_lb[j] = fluxes[j];
        p.col_ub[j] = fluxes[j];
        p.obj[j] = 0.0;
    }
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    if s.status != Status::Optimal {
        return Err(Error::Solver(format!("CycleFreeFlux ended with {:?}", s.status)));
    }
    Ok(s.x)
}
