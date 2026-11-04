//! Enzyme-constrained flux balance in the style of GECKO (Sanchez et al. 2017).
//!
//! Each enzyme `k` gets a usage variable `e_k` in mmol/gDW. A reaction's flux in a direction is bounded
//! by what its isozymes can catalyse, `|v_i| <= sum_k kcat_ik e_ik`, the total mass used is bounded by
//! the available protein pool, `sum_k MW_k e_k <= P`, and measured abundances bound single enzymes.
//! The turnover numbers and molecular weights are inputs; this module holds no biological data.

use super::fba::FluxSolution;
use crate::error::{Error, Result};
use crate::lp::{self, Sense};
use crate::model::Model;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// One isozyme's catalysis of a reaction.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Catalysis {
    /// Reaction identifier.
    pub reaction: String,
    /// Enzyme identifier.
    pub enzyme: String,
    /// Forward turnover number in 1/h.
    pub kcat_forward: f64,
    /// Reverse turnover number in 1/h; zero means no reverse catalysis.
    #[serde(default)]
    pub kcat_reverse: f64,
}

/// An enzyme.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Enzyme {
    /// Identifier.
    pub id: String,
    /// Molecular weight in g/mmol (equal to kDa).
    pub molecular_weight: f64,
    /// Measured abundance in mmol/gDW, when known.
    #[serde(default)]
    pub abundance: Option<f64>,
}

/// Enzyme data for a model.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct EnzymeData {
    /// Enzymes.
    pub enzymes: Vec<Enzyme>,
    /// Catalysis relations.
    pub catalysis: Vec<Catalysis>,
    /// Protein pool available to the enzymes with data, in g/gDW.
    pub protein_pool: f64,
}

/// Solution of an enzyme-constrained FBA.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EnzymeSolution {
    /// Fluxes and objective.
    pub flux: FluxSolution,
    /// Usage (mmol/gDW) per enzyme.
    pub usage: BTreeMap<String, f64>,
    /// Mass of protein used, g/gDW.
    pub pool_used: f64,
}

/// Runs enzyme-constrained FBA with the default backend.
pub fn enzyme_constrained_fba(model: &Model, data: &EnzymeData) -> Result<EnzymeSolution> {
    if data.protein_pool <= 0.0 {
        return Err(Error::Invalid("the protein pool must be positive".into()));
    }
    let backend = lp::default_backend();
    let mut p = model.to_lp();
    let n = model.reactions.len();
    let mut usage_col: BTreeMap<&str, usize> = BTreeMap::new();
    let mut pool_row = Vec::with_capacity(data.enzymes.len());
    for e in &data.enzymes {
        if e.molecular_weight <= 0.0 {
            return Err(Error::Invalid(format!("enzyme {} has no positive molecular weight", e.id)));
        }
        let col = p.add_col(0.0, 0.0, e.abundance.unwrap_or(f64::INFINITY), false);
        usage_col.insert(e.id.as_str(), col);
        pool_row.push((col, e.molecular_weight));
    }
    p.add_row(f64::NEG_INFINITY, data.protein_pool, pool_row.clone());
    let mut by_reaction: BTreeMap<usize, Vec<&Catalysis>> = BTreeMap::new();
    for c in &data.catalysis {
        by_reaction.entry(model.reaction_index(&c.reaction)?).or_default().push(c);
    }
    for (&j, list) in &by_reaction {
        // Isozymes share a reaction, so each contributes its own capacity to the flux bound.
        let mut forward = vec![(j, 1.0)];
        let mut reverse = vec![(j, 1.0)];
        for c in list {
            let col = *usage_col
                .get(c.enzyme.as_str())
                .ok_or_else(|| Error::Unknown { kind: "enzyme", id: c.enzyme.clone() })?;
            forward.push((col, -c.kcat_forward));
            reverse.push((col, c.kcat_reverse));
        }
        p.add_row(f64::NEG_INFINITY, 0.0, forward);
        p.add_row(0.0, f64::INFINITY, reverse);
    }
    p.sense = Some(model.sense.unwrap_or(Sense::Maximize));
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    if !s.is_optimal() {
        return Err(Error::Solver(format!("enzyme-constrained FBA ended with {:?}", s.status)));
    }
    let usage: BTreeMap<String, f64> = usage_col.iter().map(|(k, &c)| ((*k).to_string(), s.x[c])).collect();
    let pool_used = pool_row.iter().map(|&(c, w)| s.x[c] * w).sum();
    Ok(EnzymeSolution {
        flux: FluxSolution {
            status: s.status,
            objective: s.objective,
            fluxes: s.x[..n].to_vec(),
            shadow_prices: Vec::new(),
            reduced_costs: Vec::new(),
        },
        usage,
        pool_used,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Metabolite, Reaction};

    /// `-> A`, then A leaves by a fast and costly route or a slow and cheap one.
    fn model() -> Model {
        let mut m = Model::new("routes");
        m.metabolites.push(Metabolite { id: "A".into(), compartment: "c".into(), ..Metabolite::default() });
        let reaction = |id: &str, c: f64| Reaction { id: id.into(), lb: 0.0, ub: 1000.0, stoichiometry: vec![(0, c)], ..Reaction::default() };
        m.reactions.push(reaction("in", 1.0));
        m.reactions.push(reaction("route", -1.0));
        m.objective = vec![(1, 1.0)];
        m.sense = Some(Sense::Maximize);
        m.reindex();
        m
    }

    fn data(pool: f64) -> EnzymeData {
        EnzymeData {
            enzymes: vec![Enzyme { id: "E".into(), molecular_weight: 50.0, abundance: None }],
            catalysis: vec![Catalysis { reaction: "route".into(), enzyme: "E".into(), kcat_forward: 3600.0, kcat_reverse: 0.0 }],
            protein_pool: pool,
        }
    }

    #[test]
    fn pool_sets_the_flux_exactly() {
        // v = kcat e and MW e = P give v = kcat P / MW.
        let s = enzyme_constrained_fba(&model(), &data(0.1)).unwrap();
        assert!((s.flux.objective - 3600.0 * 0.1 / 50.0).abs() < 1e-8, "{}", s.flux.objective);
        assert!((s.pool_used - 0.1).abs() < 1e-9);
    }

    #[test]
    fn abundance_bound_tightens_the_limit() {
        let mut d = data(1.0);
        d.enzymes[0].abundance = Some(0.001);
        let s = enzyme_constrained_fba(&model(), &d).unwrap();
        assert!((s.flux.objective - 3.6).abs() < 1e-8);
    }

    #[test]
    fn isozymes_add_capacity() {
        let mut d = data(0.1);
        d.enzymes.push(Enzyme { id: "E2".into(), molecular_weight: 25.0, abundance: None });
        d.catalysis.push(Catalysis { reaction: "route".into(), enzyme: "E2".into(), kcat_forward: 3600.0, kcat_reverse: 0.0 });
        let s = enzyme_constrained_fba(&model(), &d).unwrap();
        // The cheaper isozyme takes all of the pool.
        assert!((s.flux.objective - 3600.0 * 0.1 / 25.0).abs() < 1e-8);
        assert!(s.usage["E"] < 1e-9);
    }
}
