//! Parsimonious FBA (Lewis et al. 2010): among fluxes keeping the objective at a fraction of its
//! optimum, minimize the total absolute flux. Each reaction is split into forward and reverse parts.

use super::fba::FluxSolution;
use crate::error::{Error, Result};
use crate::lp::{self, Backend, Problem, Sense, Status};
use crate::model::Model;

/// Runs pFBA with the default backend. `fraction` is the retained share of the optimum.
pub fn pfba(model: &Model, fraction: f64) -> Result<FluxSolution> {
    pfba_with(model, fraction, lp::default_backend().as_ref())
}

/// Runs pFBA with a chosen backend.
pub fn pfba_with(model: &Model, fraction: f64, backend: &dyn Backend) -> Result<FluxSolution> {
    let optimum = super::fba::fba_with(model, backend)?;
    if optimum.status != Status::Optimal {
        return Ok(optimum);
    }
    let n = model.reactions.len();
    let mut p = Problem { sense: Some(Sense::Minimize), ..Problem::default() };
    for r in &model.reactions {
        p.add_col(1.0, r.lb.max(0.0), r.ub.max(0.0), false);
    }
    for r in &model.reactions {
        p.add_col(1.0, (-r.ub).max(0.0), (-r.lb).max(0.0), false);
    }
    let lp_base = model.to_lp();
    for row in lp_base.rows {
        let mut coefs = Vec::with_capacity(row.len() * 2);
        for (j, c) in row {
            coefs.push((j, c));
            coefs.push((n + j, -c));
        }
        p.add_row(0.0, 0.0, coefs);
    }
    let mut coefs = Vec::new();
    for &(j, c) in &model.objective {
        coefs.push((j, c));
        coefs.push((n + j, -c));
    }
    let (lo, hi) = match model.sense.unwrap_or(Sense::Maximize) {
        Sense::Maximize => (fraction * optimum.objective, f64::INFINITY),
        Sense::Minimize => (f64::NEG_INFINITY, optimum.objective / fraction.max(1e-300)),
    };
    p.add_row(lo, hi, coefs);
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    if !s.is_optimal() {
        return Err(Error::Solver(format!("pFBA minimization ended with {:?}", s.status)));
    }
    let fluxes: Vec<f64> = (0..n).map(|j| s.x[j] - s.x[n + j]).collect();
    Ok(FluxSolution {
        status: s.status,
        objective: s.objective,
        fluxes,
        shadow_prices: Vec::new(),
        reduced_costs: Vec::new(),
    })
}
