//! Enzyme-constrained flux balance in the style of GECKO (Sanchez et al. 2017).
//!
//! Each catalysis `c` (one reaction by one isozyme) gets an allocation `a_c` in mmol/gDW. A reaction's
//! flux is bounded by what its isozymes can catalyse, `-sum kcat-_c a_c <= v <= sum kcat+_c a_c`. Each
//! enzyme is shared among the catalyses that use it, `sum_c n_ck a_c <= e_k`, with `n_ck` the copies of
//! enzyme `k` in the complex of `c`. The total mass used is bounded by the protein pool,
//! `sum_k MW_k e_k <= P`, and measured abundances bound single enzymes. The turnover numbers and
//! molecular weights are inputs; this module holds no biological data.

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
    /// Enzyme identifier, for a catalyst of one enzyme. Leave empty when `subunits` is given.
    #[serde(default)]
    pub enzyme: String,
    /// Enzymes and copies of each in a complex, all of which one catalytic cycle occupies.
    #[serde(default)]
    pub subunits: Vec<(String, f64)>,
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
    let mut by_reaction: BTreeMap<usize, Vec<(usize, &Catalysis)>> = BTreeMap::new();
    let mut shares: BTreeMap<usize, Vec<(usize, f64)>> = BTreeMap::new();
    for c in &data.catalysis {
        let j = model.reaction_index(&c.reaction)?;
        let allocation = p.add_col(0.0, 0.0, f64::INFINITY, false);
        by_reaction.entry(j).or_default().push((allocation, c));
        let single = [(c.enzyme.clone(), 1.0)];
        let members: &[(String, f64)] = if c.subunits.is_empty() { &single } else { &c.subunits };
        for (enzyme, copies) in members {
            let col = *usage_col.get(enzyme.as_str()).ok_or_else(|| Error::Unknown { kind: "enzyme", id: enzyme.clone() })?;
            shares.entry(col).or_default().push((allocation, *copies));
        }
    }
    for (&col, list) in &shares {
        let mut row = list.clone();
        row.push((col, -1.0));
        p.add_row(f64::NEG_INFINITY, 0.0, row);
    }
    for (&j, list) in &by_reaction {
        // Isozymes share a reaction, so each contributes its own capacity to the flux bound.
        let mut forward = vec![(j, 1.0)];
        let mut reverse = vec![(j, 1.0)];
        for &(allocation, c) in list {
            forward.push((allocation, -c.kcat_forward));
            reverse.push((allocation, c.kcat_reverse));
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
        let reaction = |id: &str, c: f64| Reaction {
            id: id.into(),
            lb: 0.0,
            ub: 1000.0,
            stoichiometry: vec![(0, c)],
            ..Reaction::default()
        };
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
            catalysis: vec![Catalysis {
                reaction: "route".into(),
                enzyme: "E".into(),
                subunits: Vec::new(),
                kcat_forward: 3600.0,
                kcat_reverse: 0.0,
            }],
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
        d.catalysis.push(Catalysis {
            reaction: "route".into(),
            enzyme: "E2".into(),
            subunits: Vec::new(),
            kcat_forward: 3600.0,
            kcat_reverse: 0.0,
        });
        let s = enzyme_constrained_fba(&model(), &d).unwrap();
        // The cheaper isozyme takes all of the pool.
        assert!((s.flux.objective - 3600.0 * 0.1 / 25.0).abs() < 1e-8);
        assert!(s.usage["E"] < 1e-9);
    }

    /// `-> A`, with A leaving by two parallel routes.
    fn parallel() -> Model {
        let mut m = model();
        m.reactions.push(Reaction { id: "other".into(), lb: 0.0, ub: 1000.0, stoichiometry: vec![(0, -1.0)], ..Reaction::default() });
        m.objective = vec![(1, 1.0), (2, 1.0)];
        m.reindex();
        m
    }

    fn catalysis(reaction: &str, enzyme: &str, kcat: f64) -> Catalysis {
        Catalysis { reaction: reaction.into(), enzyme: enzyme.into(), subunits: Vec::new(), kcat_forward: kcat, kcat_reverse: 0.0 }
    }

    #[test]
    fn a_promiscuous_enzyme_shares_its_capacity() {
        let d = EnzymeData {
            enzymes: vec![Enzyme { id: "E".into(), molecular_weight: 1.0, abundance: None }],
            catalysis: vec![catalysis("route", "E", 10.0), catalysis("other", "E", 10.0)],
            protein_pool: 1.0,
        };
        let s = enzyme_constrained_fba(&parallel(), &d).unwrap();
        assert!((s.flux.objective - 10.0).abs() < 1e-8, "{}", s.flux.objective);
        assert!((s.usage["E"] - 1.0).abs() < 1e-8);
    }

    #[test]
    fn a_complex_takes_all_of_its_subunits() {
        // One cycle of the complex takes one A and two B, of 10 and 5 g/mmol, so it costs 20 g/mmol.
        let d = EnzymeData {
            enzymes: vec![
                Enzyme { id: "A".into(), molecular_weight: 10.0, abundance: None },
                Enzyme { id: "B".into(), molecular_weight: 5.0, abundance: None },
            ],
            catalysis: vec![Catalysis {
                reaction: "route".into(),
                enzyme: String::new(),
                subunits: vec![("A".into(), 1.0), ("B".into(), 2.0)],
                kcat_forward: 100.0,
                kcat_reverse: 0.0,
            }],
            protein_pool: 2.0,
        };
        let s = enzyme_constrained_fba(&model(), &d).unwrap();
        assert!((s.flux.objective - 100.0 * 2.0 / 20.0).abs() < 1e-8, "{}", s.flux.objective);
        assert!((s.usage["B"] - 2.0 * s.usage["A"]).abs() < 1e-8);
    }

    #[test]
    fn abundance_of_one_subunit_limits_the_complex() {
        let d = EnzymeData {
            enzymes: vec![
                Enzyme { id: "A".into(), molecular_weight: 1.0, abundance: None },
                Enzyme { id: "B".into(), molecular_weight: 1.0, abundance: Some(0.01) },
            ],
            catalysis: vec![Catalysis {
                reaction: "route".into(),
                enzyme: String::new(),
                subunits: vec![("A".into(), 1.0), ("B".into(), 1.0)],
                kcat_forward: 100.0,
                kcat_reverse: 0.0,
            }],
            protein_pool: 10.0,
        };
        let s = enzyme_constrained_fba(&model(), &d).unwrap();
        assert!((s.flux.objective - 1.0).abs() < 1e-8, "{}", s.flux.objective);
    }
}
