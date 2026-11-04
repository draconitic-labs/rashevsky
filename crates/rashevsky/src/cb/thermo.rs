//! Thermodynamics-based flux analysis (Henry, Broadbelt and Hatzimanikatis 2007).
//!
//! Each reaction with a known standard transformed Gibbs energy gets a variable
//! `DrG' = DrG'0 + e + RT sum_j n_ij ln c_j`, with `e` within `k sigma` of the estimate's uncertainty,
//! the logarithms of concentrations bounded, and forward and reverse fluxes tied to the sign of
//! `DrG'` by binary variables. Reactions without data stay unconstrained. Standard energies come from
//! outside this crate: `scripts/export_thermo.py` reads them from eQuilibrator (component
//! contribution; Noor et al. 2013, Beber et al. 2022). Water and protons do not enter the concentration
//! term, since their activities are part of the transformed standard state.

use super::fba::FluxSolution;
use crate::error::{Error, Result};
use crate::lp::{self, Sense};
use crate::model::Model;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Gas constant in kJ/(mol K).
pub const GAS_CONSTANT: f64 = 8.314_462_618e-3;

/// Standard transformed Gibbs energy of a reaction.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ReactionEnergy {
    /// `DrG'0` in kJ/mol.
    pub dg0: f64,
    /// Standard uncertainty in kJ/mol.
    pub sigma: f64,
}

/// Thermodynamic data for a model.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ThermoData {
    /// Temperature in K.
    pub temperature: f64,
    /// Energies by reaction identifier.
    pub reactions: BTreeMap<String, ReactionEnergy>,
    /// Metabolites whose concentration term is omitted, such as water and protons.
    pub excluded_metabolites: Vec<String>,
    /// Concentration bounds (M) by metabolite identifier.
    pub concentration_bounds: BTreeMap<String, (f64, f64)>,
    /// Source of the numbers.
    pub source: String,
}

/// Settings of [`tfa`].
#[derive(Clone, Copy, Debug)]
pub struct TfaOptions {
    /// Lower concentration bound (M) for metabolites without their own.
    pub c_min: f64,
    /// Upper concentration bound (M) for metabolites without their own.
    pub c_max: f64,
    /// Allowed error in units of the estimate's standard uncertainty.
    pub sigmas: f64,
    /// Bound on `|DrG'|` in the big-M rows, kJ/mol.
    pub big_m: f64,
    /// Smallest `|DrG'|` (kJ/mol) that counts as a driving force.
    pub epsilon: f64,
}

impl Default for TfaOptions {
    fn default() -> Self {
        Self { c_min: 1e-6, c_max: 0.02, sigmas: 2.0, big_m: 1000.0, epsilon: 1e-3 }
    }
}

/// Solution of a thermodynamics-based flux analysis.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TfaSolution {
    /// Fluxes and objective.
    pub flux: FluxSolution,
    /// Concentration (M) by metabolite identifier, for metabolites in the constraints.
    pub concentrations: BTreeMap<String, f64>,
    /// `DrG'` (kJ/mol) by reaction identifier.
    pub gibbs: BTreeMap<String, f64>,
}

/// Runs TFA with the default backend.
pub fn tfa(model: &Model, data: &ThermoData, options: &TfaOptions) -> Result<TfaSolution> {
    let backend = lp::default_backend();
    let rt = GAS_CONSTANT * if data.temperature > 0.0 { data.temperature } else { 298.15 };
    let mut p = model.to_lp();
    let n = model.reactions.len();
    let bound = |id: &str| -> (f64, f64) {
        data.concentration_bounds.get(id).copied().unwrap_or((options.c_min, options.c_max))
    };
    let bounds = super::finite_bounds(model, lp::default_backend().as_ref())?;
    // Logarithm of concentration per metabolite that appears in a constrained reaction.
    let mut log_col: BTreeMap<usize, usize> = BTreeMap::new();
    let mut gibbs_col: BTreeMap<usize, usize> = BTreeMap::new();
    for (id, energy) in &data.reactions {
        let j = model.reaction_index(id)?;
        let r = &model.reactions[j];
        let mut coefs: Vec<(usize, f64)> = Vec::new();
        for &(i, nu) in &r.stoichiometry {
            let met = &model.metabolites[i];
            if data.excluded_metabolites.iter().any(|e| e == &met.id) {
                continue;
            }
            let col = *log_col.entry(i).or_insert_with(|| {
                let (lo, hi) = bound(&met.id);
                p.add_col(0.0, lo.ln(), hi.ln(), false)
            });
            coefs.push((col, nu * rt));
        }
        let m = options.big_m;
        let g = p.add_col(0.0, -m, m, false);
        gibbs_col.insert(j, g);
        let error = p.add_col(0.0, -options.sigmas * energy.sigma, options.sigmas * energy.sigma, false);
        // DrG' - error - RT sum n ln c = DrG'0
        coefs.push((g, -1.0));
        coefs.push((error, 1.0));
        p.add_row(-energy.dg0, -energy.dg0, coefs);
        // Forward flux needs DrG' <= -eps, reverse flux needs DrG' >= eps.
        let (lb, ub) = bounds[j];
        let zf = p.add_col(0.0, 0.0, 1.0, true);
        let zr = p.add_col(0.0, 0.0, 1.0, true);
        p.add_row(f64::NEG_INFINITY, 0.0, vec![(j, 1.0), (zf, -ub.max(0.0))]);
        p.add_row(0.0, f64::INFINITY, vec![(j, 1.0), (zr, -lb.min(0.0))]);
        p.add_row(f64::NEG_INFINITY, m - options.epsilon, vec![(g, 1.0), (zf, m)]);
        p.add_row(options.epsilon - m, f64::INFINITY, vec![(g, 1.0), (zr, -m)]);
        p.add_row(f64::NEG_INFINITY, 1.0, vec![(zf, 1.0), (zr, 1.0)]);
    }
    p.sense = Some(model.sense.unwrap_or(Sense::Maximize));
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    if !s.is_optimal() {
        return Err(Error::Solver(format!("TFA ended with {:?}", s.status)));
    }
    let concentrations = log_col.iter().map(|(&i, &c)| (model.metabolites[i].id.clone(), s.x[c].exp())).collect();
    let gibbs = gibbs_col.iter().map(|(&j, &g)| (model.reactions[j].id.clone(), s.x[g])).collect();
    Ok(TfaSolution {
        flux: FluxSolution {
            status: s.status,
            objective: s.objective,
            fluxes: s.x[..n].to_vec(),
            shadow_prices: Vec::new(),
            reduced_costs: Vec::new(),
        },
        concentrations,
        gibbs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Metabolite, Reaction};

    /// `-> A`, `A -> B`, `B ->`, with concentrations bounded so the middle step is fixed in sign.
    fn chain() -> Model {
        let mut m = Model::new("chain");
        for id in ["A", "B"] {
            m.metabolites.push(Metabolite { id: id.into(), compartment: "c".into(), ..Metabolite::default() });
        }
        let reaction = |id: &str, st: Vec<(usize, f64)>| Reaction {
            id: id.into(),
            lb: -10.0,
            ub: 10.0,
            stoichiometry: st,
            ..Reaction::default()
        };
        m.reactions.push(reaction("in", vec![(0, 1.0)]));
        m.reactions.push(reaction("step", vec![(0, -1.0), (1, 1.0)]));
        m.reactions.push(reaction("out", vec![(1, -1.0)]));
        m.objective = vec![(2, 1.0)];
        m.sense = Some(Sense::Maximize);
        m.reindex();
        m
    }

    fn data(dg0: f64, sigma: f64) -> ThermoData {
        let mut d = ThermoData { temperature: 298.15, ..ThermoData::default() };
        d.reactions.insert("step".into(), ReactionEnergy { dg0, sigma });
        d
    }

    #[test]
    fn favourable_reaction_keeps_its_flux() {
        let s = tfa(&chain(), &data(-30.0, 0.0), &TfaOptions::default()).unwrap();
        assert!((s.flux.objective - 10.0).abs() < 1e-6);
        assert!(s.gibbs["step"] < 0.0);
    }

    #[test]
    fn uphill_reaction_blocks_forward_flux() {
        // RT ln(0.02 / 1e-6) is 24.5 kJ/mol at 298 K, less than the +30 required.
        let s = tfa(&chain(), &data(30.0, 0.0), &TfaOptions::default()).unwrap();
        assert!(s.flux.objective <= 1e-6, "{}", s.flux.objective);
    }

    #[test]
    fn concentrations_can_overcome_a_small_barrier() {
        // +20 kJ/mol is reachable, and the solution must use a concentration ratio that does it.
        let s = tfa(&chain(), &data(20.0, 0.0), &TfaOptions::default()).unwrap();
        assert!((s.flux.objective - 10.0).abs() < 1e-6);
        let ratio = s.concentrations["B"] / s.concentrations["A"];
        let rt = GAS_CONSTANT * 298.15;
        assert!(20.0 + rt * ratio.ln() < 0.0, "ratio {ratio}");
    }

    #[test]
    fn uncertainty_widens_what_is_feasible() {
        // Within two sigma of +30 lies +20, which is reachable.
        let s = tfa(&chain(), &data(30.0, 5.0), &TfaOptions::default()).unwrap();
        assert!((s.flux.objective - 10.0).abs() < 1e-6);
    }
}
