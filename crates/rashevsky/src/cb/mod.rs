//! Constraint-based analyses over a stoichiometric [`Model`](crate::model::Model).

pub mod deletion;
pub mod dfba;
pub mod envelope;
pub mod enzyme;
pub mod fba;
pub mod fva;
pub mod loopless;
pub mod media;
pub mod moma;
pub mod pfba;
pub mod sampling;
pub mod thermo;

use crate::error::{Error, Result};
use crate::lp::{Backend, Sense, Status};
use crate::model::Model;

/// The range the objective may take when `fraction` of its optimum is retained. A fraction of zero
/// retains nothing. The limit is `fraction * optimum` or `optimum / fraction`, whichever lets the
/// objective worsen, so that it also holds for a negative optimum.
pub(crate) fn retained_range(sense: Sense, optimum: f64, fraction: f64) -> Option<(f64, f64)> {
    if fraction <= 0.0 {
        return None;
    }
    let (a, b) = (fraction * optimum, optimum / fraction);
    Some(match sense {
        Sense::Maximize => (a.min(b), f64::INFINITY),
        Sense::Minimize => (f64::NEG_INFINITY, a.max(b)),
    })
}

/// Flux bounds with each infinite side replaced by the most the network allows, for mixed-integer rows
/// that multiply a bound by a binary variable.
pub(crate) fn finite_bounds(model: &Model, backend: &dyn Backend) -> Result<Vec<(f64, f64)>> {
    let mut out: Vec<(f64, f64)> = model.reactions.iter().map(|r| (r.lb, r.ub)).collect();
    if out.iter().all(|(lb, ub)| lb.is_finite() && ub.is_finite()) {
        return Ok(out);
    }
    let mut problem = model.to_lp();
    problem.obj.iter_mut().for_each(|c| *c = 0.0);
    let mut session = backend.open(&problem)?;
    let mut obj = vec![0.0; model.reactions.len()];
    for j in 0..out.len() {
        for (infinite, sense) in [(!out[j].0.is_finite(), Sense::Minimize), (!out[j].1.is_finite(), Sense::Maximize)] {
            if !infinite {
                continue;
            }
            obj[j] = 1.0;
            session.set_objective(sense, &obj)?;
            obj[j] = 0.0;
            let s = session.solve()?;
            if s.status != Status::Optimal {
                return Err(Error::Invalid(format!(
                    "reaction {} has no finite limit ({:?}); this analysis needs one, so give it a bound",
                    model.reactions[j].id, s.status
                )));
            }
            if sense == Sense::Minimize {
                out[j].0 = s.objective;
            } else {
                out[j].1 = s.objective;
            }
        }
    }
    Ok(out)
}

pub use deletion::{
    Deletion, double_gene_deletion, double_reaction_deletion, single_gene_deletion, single_gene_deletion_with,
    single_reaction_deletion, synthetic_lethal_genes,
};
pub use dfba::{DfbaOptions, DfbaResult, Substrate, dfba};
pub use envelope::{EnvelopePoint, production_envelope};
pub use enzyme::{EnzymeData, EnzymeSolution, enzyme_constrained_fba};
pub use fba::{FluxSolution, fba, fba_with};
pub use fva::{FvaOptions, FvaRange, fva, fva_with};
pub use loopless::{cycle_free, loopless_fba};
pub use media::{MinimalMedium, medium, minimal_medium, set_medium};
pub use moma::{RoomOptions, linear_moma, moma, room};
pub use pfba::pfba;
pub use sampling::{FluxSamples, SamplerOptions, chrr, optgp};
pub use thermo::{TfaOptions, TfaSolution, ThermoData, tfa};
