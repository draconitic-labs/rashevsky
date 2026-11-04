//! Production envelopes: the range of one flux over a grid of another (the phenotypic phase plane).

use crate::error::{Error, Result};
use crate::lp::{self, Sense, Status};
use crate::model::Model;
use serde::{Deserialize, Serialize};

/// One point of a production envelope.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EnvelopePoint {
    /// Value of the controlled flux.
    pub x: f64,
    /// Smallest feasible value of the target flux, NaN if infeasible.
    pub minimum: f64,
    /// Largest feasible value of the target flux, NaN if infeasible.
    pub maximum: f64,
}

/// Computes the envelope of `target` against `control` on an evenly spaced grid of `points`.
///
/// The grid spans the feasible range of `control` unless `range` overrides it.
pub fn production_envelope(
    model: &Model,
    control: &str,
    target: &str,
    points: usize,
    range: Option<(f64, f64)>,
) -> Result<Vec<EnvelopePoint>> {
    if points < 2 {
        return Err(Error::Invalid("an envelope needs at least two points".into()));
    }
    let cj = model.reaction_index(control)?;
    let tj = model.reaction_index(target)?;
    let backend = lp::default_backend();
    let mut problem = model.to_lp();
    problem.obj.iter_mut().for_each(|c| *c = 0.0);
    let mut session = backend.open(&problem)?;
    let mut obj = vec![0.0; model.reactions.len()];
    let mut extreme = |session: &mut Box<dyn lp::Session>, j: usize, sense: Sense| -> Result<Option<f64>> {
        obj[j] = 1.0;
        session.set_objective(sense, &obj)?;
        obj[j] = 0.0;
        let s = session.solve()?;
        Ok(match s.status {
            Status::Optimal => Some(s.objective),
            _ => None,
        })
    };
    let (lo, hi) = match range {
        Some(r) => r,
        None => {
            let lo = extreme(&mut session, cj, Sense::Minimize)?;
            let hi = extreme(&mut session, cj, Sense::Maximize)?;
            match (lo, hi) {
                (Some(a), Some(b)) => (a, b),
                _ => return Err(Error::Solver("the control flux has no feasible range".into())),
            }
        }
    };
    let mut out = Vec::with_capacity(points);
    for k in 0..points {
        let x = lo + (hi - lo) * k as f64 / (points - 1) as f64;
        session.set_col_bounds(cj, x, x)?;
        let minimum = extreme(&mut session, tj, Sense::Minimize)?.unwrap_or(f64::NAN);
        let maximum = extreme(&mut session, tj, Sense::Maximize)?.unwrap_or(f64::NAN);
        out.push(EnvelopePoint { x, minimum, maximum });
    }
    Ok(out)
}
