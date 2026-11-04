//! Flux balance analysis.
//!
//! `max c'v` subject to `S v = 0` and `lb <= v <= ub` (Orth, Thiele and Palsson 2010).

use crate::error::{Error, Result};
use crate::lp::{self, Backend, Status};
use crate::model::Model;
use serde::{Deserialize, Serialize};

/// A flux distribution with its duals.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FluxSolution {
    /// Solver status.
    pub status: Status,
    /// Objective value, NaN unless optimal.
    pub objective: f64,
    /// Flux per reaction.
    pub fluxes: Vec<f64>,
    /// Shadow price per metabolite.
    pub shadow_prices: Vec<f64>,
    /// Reduced cost per reaction.
    pub reduced_costs: Vec<f64>,
}

impl FluxSolution {
    /// Flux of a reaction by identifier.
    pub fn flux(&self, model: &Model, id: &str) -> Result<f64> {
        Ok(self.fluxes[model.reaction_index(id)?])
    }
}

/// Runs FBA with the default backend.
pub fn fba(model: &Model) -> Result<FluxSolution> {
    fba_with(model, lp::default_backend().as_ref())
}

/// Runs FBA with a chosen backend.
pub fn fba_with(model: &Model, backend: &dyn Backend) -> Result<FluxSolution> {
    if model.objective.is_empty() {
        return Err(Error::Invalid("the model has no objective".into()));
    }
    let mut session = backend.open(&model.to_lp())?;
    let s = session.solve()?;
    Ok(FluxSolution {
        status: s.status,
        objective: s.objective,
        fluxes: s.x,
        shadow_prices: s.row_dual,
        reduced_costs: s.col_dual,
    })
}
