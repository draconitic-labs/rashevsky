//! Flux variability analysis (Mahadevan and Schilling 2003).
//!
//! The objective is first optimized, then held at a fraction of its optimum by one extra row, and
//! each reaction is minimized and maximized over the remaining feasible set. Each worker keeps one
//! solver session, so every LP restarts from the previous basis (Gudmundsson and Thiele 2010).

use crate::error::{Error, Result};
use crate::lp::{self, Backend, Sense, Session};
use crate::model::Model;
use serde::{Deserialize, Serialize};

/// Settings of [`fva`].
#[derive(Clone, Debug)]
pub struct FvaOptions {
    /// Fraction of the optimal objective that must be retained, in `[0, 1]`.
    pub fraction_of_optimum: f64,
    /// Reactions to analyse; `None` means all of them.
    pub reactions: Option<Vec<usize>>,
}

impl Default for FvaOptions {
    fn default() -> Self {
        Self { fraction_of_optimum: 1.0, reactions: None }
    }
}

/// Flux range of one reaction.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FvaRange {
    /// Reaction index.
    pub reaction: usize,
    /// Reaction identifier.
    pub id: String,
    /// Minimum flux.
    pub minimum: f64,
    /// Maximum flux.
    pub maximum: f64,
}

/// Extends the steady-state LP of `model` with the objective-retention row.
/// Returns the problem with its objective cleared, ready for per-reaction objectives.
pub(crate) fn retained_problem(
    model: &Model,
    backend: &dyn Backend,
    fraction: f64,
) -> Result<(lp::Problem, f64)> {
    if !(0.0..=1.0).contains(&fraction) {
        return Err(Error::Invalid("fraction_of_optimum must lie in [0, 1]".into()));
    }
    let mut problem = model.to_lp();
    let optimum = {
        let mut session = backend.open(&problem)?;
        let s = session.solve()?;
        if !s.is_optimal() {
            return Err(Error::Solver(format!("FBA status {:?}", s.status)));
        }
        s.objective
    };
    let coefs: Vec<(usize, f64)> = problem
        .obj
        .iter()
        .enumerate()
        .filter(|(_, c)| **c != 0.0)
        .map(|(j, c)| (j, *c))
        .collect();
    // Loosen by the solver's feasibility tolerance so that a fraction of 1 stays feasible.
    let margin = 1e-9 * (1.0 + optimum.abs());
    let (lb, ub) = match model.sense.unwrap_or(Sense::Maximize) {
        Sense::Maximize => (fraction * optimum - margin, f64::INFINITY),
        Sense::Minimize => (f64::NEG_INFINITY, optimum / fraction.max(1e-300) + margin),
    };
    problem.add_row(lb, ub, coefs);
    problem.obj.iter_mut().for_each(|c| *c = 0.0);
    Ok((problem, optimum))
}

fn range_of(session: &mut dyn Session, j: usize, obj: &mut [f64]) -> Result<(f64, f64)> {
    obj[j] = 1.0;
    let mut out = [0.0; 2];
    for (k, sense) in [Sense::Minimize, Sense::Maximize].into_iter().enumerate() {
        session.set_objective(sense, obj)?;
        let s = session.solve()?;
        out[k] = match s.status {
            lp::Status::Optimal => s.objective,
            lp::Status::Unbounded => {
                if sense == Sense::Minimize { f64::NEG_INFINITY } else { f64::INFINITY }
            }
            other => return Err(Error::Solver(format!("FVA LP on column {j} ended with {other:?}"))),
        };
    }
    obj[j] = 0.0;
    Ok((out[0], out[1]))
}

/// Flux variability analysis with the default backend.
pub fn fva(model: &Model, options: &FvaOptions) -> Result<Vec<FvaRange>> {
    fva_with(model, options, lp::default_backend().as_ref())
}

/// Flux variability analysis with a chosen backend.
pub fn fva_with(model: &Model, options: &FvaOptions, backend: &dyn Backend) -> Result<Vec<FvaRange>> {
    let (problem, _) = retained_problem(model, backend, options.fraction_of_optimum)?;
    let targets: Vec<usize> =
        options.reactions.clone().unwrap_or_else(|| (0..model.reactions.len()).collect());
    let chunk = crate::parallel::chunk_size(targets.len());
    let parts: Vec<Result<Vec<FvaRange>>> = crate::parallel::map_chunks(&targets, chunk, |slice| {
        let mut session = backend.open(&problem)?;
        let mut obj = vec![0.0; problem.n_cols()];
        slice
            .iter()
            .map(|&j| {
                let (minimum, maximum) = range_of(session.as_mut(), j, &mut obj)?;
                Ok(FvaRange { reaction: j, id: model.reactions[j].id.clone(), minimum, maximum })
            })
            .collect()
    });
    let mut out = Vec::with_capacity(targets.len());
    for part in parts {
        out.extend(part?);
    }
    Ok(out)
}
