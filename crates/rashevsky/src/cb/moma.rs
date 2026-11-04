//! Responses of a perturbed network relative to a reference flux distribution.
//!
//! MOMA (Segre, Vitkup and Church 2002) minimizes the Euclidean distance to the reference; the
//! linear variant minimizes the L1 distance. ROOM (Shlomi, Berkman and Ruppin 2005) minimizes the
//! number of fluxes that leave a tolerance band around the reference.

use super::fba::FluxSolution;
use crate::error::{Error, Result};
use crate::lp::{self, Backend, Problem, Sense};
use crate::model::Model;

/// MOMA: minimize `sum (v - w)^2` with a quadratic program.
pub fn moma(model: &Model, reference: &[f64]) -> Result<FluxSolution> {
    moma_with(model, reference, lp::default_backend().as_ref())
}

/// MOMA with a chosen backend, which must solve convex quadratic programs.
pub fn moma_with(model: &Model, reference: &[f64], backend: &dyn Backend) -> Result<FluxSolution> {
    check_len(model, reference)?;
    let mut p = model.to_lp();
    p.sense = Some(Sense::Minimize);
    let n = model.reactions.len();
    for j in 0..n {
        p.obj[j] = -2.0 * reference[j];
    }
    p.hessian = (0..n).map(|j| (j, j, 2.0)).collect();
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    let distance: f64 = if s.is_optimal() {
        s.x.iter().zip(reference).map(|(v, w)| (v - w) * (v - w)).sum()
    } else {
        f64::NAN
    };
    Ok(FluxSolution {
        status: s.status,
        objective: distance,
        fluxes: s.x,
        shadow_prices: Vec::new(),
        reduced_costs: Vec::new(),
    })
}

/// Linear MOMA: minimize `sum |v - w|`.
pub fn linear_moma(model: &Model, reference: &[f64]) -> Result<FluxSolution> {
    linear_moma_with(model, reference, lp::default_backend().as_ref())
}

/// Linear MOMA with a chosen backend.
pub fn linear_moma_with(model: &Model, reference: &[f64], backend: &dyn Backend) -> Result<FluxSolution> {
    check_len(model, reference)?;
    let mut p = model.to_lp();
    p.sense = Some(Sense::Minimize);
    let n = model.reactions.len();
    p.obj.iter_mut().for_each(|c| *c = 0.0);
    for j in 0..n {
        let d = p.add_col(1.0, 0.0, f64::INFINITY, false);
        p.add_row(f64::NEG_INFINITY, reference[j], vec![(j, 1.0), (d, -1.0)]);
        p.add_row(reference[j], f64::INFINITY, vec![(j, 1.0), (d, 1.0)]);
    }
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    Ok(FluxSolution {
        status: s.status,
        objective: s.objective,
        fluxes: s.x[..n].to_vec(),
        shadow_prices: Vec::new(),
        reduced_costs: Vec::new(),
    })
}

/// Settings of [`room`].
#[derive(Clone, Copy, Debug)]
pub struct RoomOptions {
    /// Relative width of the tolerance band.
    pub delta: f64,
    /// Absolute width of the tolerance band.
    pub epsilon: f64,
    /// Solve the linear relaxation instead of the mixed-integer program.
    pub relax: bool,
}

impl Default for RoomOptions {
    fn default() -> Self {
        Self { delta: 0.03, epsilon: 1e-3, relax: false }
    }
}

/// ROOM: minimize the number of significant flux changes.
pub fn room(model: &Model, reference: &[f64], options: RoomOptions) -> Result<FluxSolution> {
    room_with(model, reference, options, lp::default_backend().as_ref())
}

/// ROOM with a chosen backend. The mixed-integer form needs a backend that solves such programs.
pub fn room_with(
    model: &Model,
    reference: &[f64],
    options: RoomOptions,
    backend: &dyn Backend,
) -> Result<FluxSolution> {
    check_len(model, reference)?;
    let mut p: Problem = model.to_lp();
    p.sense = Some(Sense::Minimize);
    let n = model.reactions.len();
    p.obj.iter_mut().for_each(|c| *c = 0.0);
    for j in 0..n {
        let w = reference[j];
        let wu = w + options.delta * w.abs() + options.epsilon;
        let wl = w - options.delta * w.abs() - options.epsilon;
        let y = p.add_col(1.0, 0.0, 1.0, !options.relax);
        let (lb, ub) = (model.reactions[j].lb, model.reactions[j].ub);
        // y = 0 keeps v inside the band; y = 1 releases it to the reaction's bounds.
        p.add_row(f64::NEG_INFINITY, wu, vec![(j, 1.0), (y, -(ub - wu).max(0.0))]);
        p.add_row(wl, f64::INFINITY, vec![(j, 1.0), (y, (wl - lb).max(0.0))]);
    }
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    Ok(FluxSolution {
        status: s.status,
        objective: s.objective,
        fluxes: s.x[..n].to_vec(),
        shadow_prices: Vec::new(),
        reduced_costs: Vec::new(),
    })
}

fn check_len(model: &Model, reference: &[f64]) -> Result<()> {
    if reference.len() != model.reactions.len() {
        return Err(Error::Invalid("reference flux length differs from reaction count".into()));
    }
    Ok(())
}
