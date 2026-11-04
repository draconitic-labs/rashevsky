//! Constraint-based model: metabolites, reactions, genes, bounds and objective.

use crate::error::{Error, Result};
use crate::gpr::{GeneIndex, Gpr};
use crate::lp::{Problem, Sense};
use crate::sbml::{self, Document};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

/// A metabolite of a stoichiometric model.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Metabolite {
    /// Identifier with the SBML `M_` prefix removed.
    pub id: String,
    /// Name.
    pub name: String,
    /// Compartment identifier.
    pub compartment: String,
    /// Chemical formula, when known.
    pub formula: Option<String>,
    /// Formal charge, when known.
    pub charge: Option<i32>,
    /// True if the metabolite is exempt from the steady-state balance.
    pub boundary: bool,
    /// `bqbiol:is` annotation resources.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub annotations: Vec<String>,
}

/// A reaction of a stoichiometric model.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Reaction {
    /// Identifier with the SBML `R_` prefix removed.
    pub id: String,
    /// Name.
    pub name: String,
    /// Lower flux bound.
    pub lb: f64,
    /// Upper flux bound.
    pub ub: f64,
    /// `(metabolite index, coefficient)`; negative for substrates.
    pub stoichiometry: Vec<(usize, f64)>,
    /// Gene-protein-reaction rule.
    pub gpr: Gpr,
    /// Subsystem from the groups package, when present.
    pub subsystem: Option<String>,
    /// `bqbiol:is` annotation resources.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub annotations: Vec<String>,
}

impl Reaction {
    /// Whether this is a biomass pseudo-reaction, going by its name. It touches every precursor, so
    /// automatic maps leave it out.
    pub fn is_biomass(&self) -> bool {
        self.id.get(..7).is_some_and(|start| start.eq_ignore_ascii_case("biomass"))
    }
}

/// A gene.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Gene {
    /// External identifier (`fbc:label` when present).
    pub id: String,
    /// Name.
    pub name: String,
}

/// A genome-scale or smaller constraint-based model.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Model {
    /// Identifier.
    pub id: String,
    /// Name.
    pub name: String,
    /// Metabolites.
    pub metabolites: Vec<Metabolite>,
    /// Reactions.
    pub reactions: Vec<Reaction>,
    /// Genes.
    pub genes: Vec<Gene>,
    /// `(reaction index, coefficient)` of the active objective.
    pub objective: Vec<(usize, f64)>,
    /// Direction of the objective.
    pub sense: Option<Sense>,
    /// Compartment identifier to name.
    pub compartments: BTreeMap<String, String>,
    #[serde(skip)]
    rxn_index: HashMap<String, usize>,
    #[serde(skip)]
    met_index: HashMap<String, usize>,
    #[serde(skip)]
    gene_index: HashMap<String, usize>,
}

fn strip(id: &str, prefix: &str) -> String {
    id.strip_prefix(prefix).unwrap_or(id).to_string()
}

impl Model {
    /// An empty model with the given identifier.
    pub fn new(id: &str) -> Self {
        Model { id: id.to_string(), ..Model::default() }
    }

    /// Reads an SBML file into a constraint-based model.
    pub fn read_sbml(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_document(&sbml::read::read_file(path)?)
    }

    /// Reads SBML text into a constraint-based model.
    pub fn from_sbml_str(xml: &str) -> Result<Self> {
        Self::from_document(&sbml::read::read_str(xml)?)
    }

    /// Builds the model from a parsed SBML document.
    pub fn from_document(doc: &Document) -> Result<Self> {
        let m = &doc.model;
        if let Some(what) = &m.unsupported_fbc {
            return Err(Error::Invalid(format!(
                "this model uses {what} from FBC version 3, which the flux analyses do not implement"
            )));
        }
        // Initial assignments and rules can define flux bounds and stoichiometries. The kinetic
        // compiler already evaluates them, so its initial values serve when any are present.
        let has_math_stoichiometry = m
            .reactions
            .iter()
            .flat_map(|r| r.reactants.iter().chain(&r.products))
            .any(|sr| sr.stoichiometry_math.is_some());
        let evaluated = if m.initial_assignments.is_empty() && m.rules.is_empty() && !has_math_stoichiometry {
            None
        } else {
            crate::kinetic::KineticModel::compile(doc).ok()
        };
        let initial_value = |id: &str| -> Option<f64> {
            let km = evaluated.as_ref()?;
            let v = km.initial[*km.index.get(id)?];
            v.is_finite().then_some(v).or(if v.is_infinite() { Some(v) } else { None })
        };
        let params: HashMap<&str, f64> =
            m.parameters.iter().filter_map(|p| initial_value(&p.id).or(p.value).map(|v| (p.id.as_str(), v))).collect();
        let mut model = Model { id: m.id.clone(), name: m.name.clone(), ..Model::default() };
        for c in &m.compartments {
            model.compartments.insert(c.id.clone(), if c.name.is_empty() { c.id.clone() } else { c.name.clone() });
        }
        let mut species_index = HashMap::with_capacity(m.species.len());
        for (i, s) in m.species.iter().enumerate() {
            species_index.insert(s.id.as_str(), i);
            model.metabolites.push(Metabolite {
                id: strip(&s.id, "M_"),
                name: s.name.clone(),
                compartment: s.compartment.clone(),
                formula: s.formula.clone(),
                charge: s.charge,
                boundary: s.boundary_condition,
                annotations: s.annotations.clone(),
            });
        }
        // Gene identifiers: SBML id mapped to the external label.
        let mut gene_labels: HashMap<&str, String> = HashMap::new();
        for g in &m.gene_products {
            // iCH360 writes `G_b4025` in the label as well, and the prefix is not part of the name.
            let label = strip(if g.label.is_empty() { &g.id } else { &g.label }, "G_");
            gene_labels.insert(g.id.as_str(), label);
        }
        let mut genes = GeneIndex::default();
        for g in &m.gene_products {
            genes.intern(&gene_labels[g.id.as_str()]);
        }
        // Subsystems from groups: the kind `partonomy` or `collection` of reaction members.
        let mut subsystem: HashMap<&str, &str> = HashMap::new();
        for g in &m.groups {
            for member in &g.members {
                let label = if g.name.is_empty() { g.id.as_str() } else { g.name.as_str() };
                subsystem.entry(member.as_str()).or_insert(label);
            }
        }
        for r in &m.reactions {
            let mut parsed = read_reaction(r, &params, &species_index, &gene_labels, &mut genes, &subsystem)?;
            // Stoichiometries assigned by name replace the attribute's value.
            if evaluated.is_some() {
                let mut merged: BTreeMap<usize, f64> = BTreeMap::new();
                for (refs, sign) in [(&r.reactants, -1.0), (&r.products, 1.0)] {
                    for sr in refs {
                        let i = species_index[sr.species.as_str()];
                        let from_math = sr.stoichiometry_math.as_ref().and_then(|e| {
                            let km = evaluated.as_ref()?;
                            let c = crate::kinetic::compile::compile_expr(e, &km.index).ok()?;
                            Some(crate::kinetic::compile::eval(&c, &km.initial, 0.0, &crate::kinetic::compile::NoDelay))
                        });
                        let value = from_math
                            .or_else(|| sr.id.as_deref().and_then(&initial_value))
                            .or(sr.stoichiometry)
                            .unwrap_or(1.0);
                        *merged.entry(i).or_insert(0.0) += sign * value;
                    }
                }
                parsed.stoichiometry = merged.into_iter().filter(|(_, c)| *c != 0.0).collect();
            }
            model.reactions.push(parsed);
        }
        let names: HashMap<&str, &str> =
            m.gene_products.iter().map(|g| (gene_labels[g.id.as_str()].as_str(), g.name.as_str())).collect();
        model.genes = genes
            .ids
            .iter()
            .map(|id| Gene { id: id.clone(), name: names.get(id.as_str()).unwrap_or(&"").to_string() })
            .collect();
        let active = m
            .active_objective
            .as_deref()
            .and_then(|a| m.objectives.iter().find(|o| o.id == a))
            .or_else(|| m.objectives.first());
        if active.is_none() {
            // The COBRA convention before FBC: an objective coefficient among the local parameters.
            for (j, r) in m.reactions.iter().enumerate() {
                let coef = local_value(r, "OBJECTIVE_COEFFICIENT").unwrap_or(0.0);
                if coef != 0.0 {
                    model.objective.push((j, coef));
                    model.sense = Some(Sense::Maximize);
                }
            }
        }
        if let Some(obj) = active {
            model.sense = Some(if obj.maximize { Sense::Maximize } else { Sense::Minimize });
            let by_id: HashMap<&str, usize> = m.reactions.iter().enumerate().map(|(i, r)| (r.id.as_str(), i)).collect();
            for (rid, coef) in &obj.flux_objectives {
                let j = *by_id.get(rid.as_str()).ok_or_else(|| Error::Unknown { kind: "reaction", id: rid.clone() })?;
                model.objective.push((j, *coef));
            }
        }
        model.reindex();
        Ok(model)
    }

    /// Rebuilds the identifier lookup tables after direct edits of the object lists.
    pub fn reindex(&mut self) {
        self.rxn_index = self.reactions.iter().enumerate().map(|(i, r)| (r.id.clone(), i)).collect();
        self.met_index = self.metabolites.iter().enumerate().map(|(i, m)| (m.id.clone(), i)).collect();
        self.gene_index = self.genes.iter().enumerate().map(|(i, g)| (g.id.clone(), i)).collect();
    }

    /// Index of a reaction by identifier.
    pub fn reaction_index(&self, id: &str) -> Result<usize> {
        self.rxn_index.get(id).copied().ok_or_else(|| Error::Unknown { kind: "reaction", id: id.into() })
    }

    /// Index of a metabolite by identifier.
    pub fn metabolite_index(&self, id: &str) -> Result<usize> {
        self.met_index.get(id).copied().ok_or_else(|| Error::Unknown { kind: "metabolite", id: id.into() })
    }

    /// Index of a gene by identifier.
    pub fn gene_index(&self, id: &str) -> Result<usize> {
        self.gene_index.get(id).copied().ok_or_else(|| Error::Unknown { kind: "gene", id: id.into() })
    }

    /// Sets the bounds of a reaction by identifier.
    pub fn set_bounds(&mut self, id: &str, lb: f64, ub: f64) -> Result<()> {
        let j = self.reaction_index(id)?;
        self.reactions[j].lb = lb;
        self.reactions[j].ub = ub;
        Ok(())
    }

    /// Replaces the objective with a single reaction, maximized by default.
    pub fn set_objective(&mut self, id: &str, sense: Sense) -> Result<()> {
        let j = self.reaction_index(id)?;
        self.objective = vec![(j, 1.0)];
        self.sense = Some(sense);
        Ok(())
    }

    /// Objective coefficient vector over reactions.
    pub fn objective_vector(&self) -> Vec<f64> {
        let mut c = vec![0.0; self.reactions.len()];
        for &(j, coef) in &self.objective {
            c[j] += coef;
        }
        c
    }

    /// The steady-state flux linear program: `S v = 0`, `lb <= v <= ub`, optimize the objective.
    /// Rows follow the metabolite order; boundary metabolites keep an empty, always-satisfied row.
    pub fn to_lp(&self) -> Problem {
        let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); self.metabolites.len()];
        for (j, r) in self.reactions.iter().enumerate() {
            for &(i, coef) in &r.stoichiometry {
                if !self.metabolites[i].boundary {
                    rows[i].push((j, coef));
                }
            }
        }
        let m = rows.len();
        Problem {
            sense: Some(self.sense.unwrap_or(Sense::Maximize)),
            obj: self.objective_vector(),
            col_lb: self.reactions.iter().map(|r| r.lb).collect(),
            col_ub: self.reactions.iter().map(|r| r.ub).collect(),
            row_lb: vec![0.0; m],
            row_ub: vec![0.0; m],
            rows,
            integer: Vec::new(),
            hessian: Vec::new(),
        }
    }

    /// Reactions disabled when the genes in `deleted` are removed.
    pub fn reactions_disabled_by(&self, deleted: &[bool]) -> Vec<usize> {
        self.reactions
            .iter()
            .enumerate()
            .filter(|(_, r)| r.gpr != Gpr::None && !r.gpr.active(deleted))
            .map(|(j, _)| j)
            .collect()
    }

    /// Coefficients of a reaction on metabolites that are balanced, leaving out boundary species.
    fn balanced_terms(&self, j: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        self.reactions[j].stoichiometry.iter().copied().filter(|&(i, _)| !self.metabolites[i].boundary)
    }

    /// True for a reaction that crosses the system's boundary: at most one balanced metabolite, whether
    /// written `A <->` or through a boundary species as `A_b <-> A`.
    pub fn is_boundary_reaction(&self, j: usize) -> bool {
        self.balanced_terms(j).count() <= 1
    }

    /// The coefficient of an exchange reaction on its one balanced metabolite. Positive flux takes
    /// metabolite up when this is positive.
    pub fn exchange_coefficient(&self, j: usize) -> Option<f64> {
        let mut terms = self.balanced_terms(j);
        match (terms.next(), terms.next()) {
            (Some((_, c)), None) => Some(c),
            _ => None,
        }
    }

    /// Exchange reactions: those with exactly one balanced metabolite.
    pub fn exchanges(&self) -> Vec<usize> {
        (0..self.reactions.len()).filter(|&j| self.exchange_coefficient(j).is_some()).collect()
    }

    /// Size summary.
    pub fn summary(&self) -> ModelSummary {
        ModelSummary {
            id: self.id.clone(),
            name: self.name.clone(),
            reactions: self.reactions.len(),
            metabolites: self.metabolites.len(),
            genes: self.genes.len(),
            compartments: self.compartments.len(),
            exchanges: self.exchanges().len(),
            objective: self.objective.iter().map(|&(j, c)| (self.reactions[j].id.clone(), c)).collect(),
        }
    }

    /// Element and charge balance of one reaction. Empty means balanced or undetermined by data.
    pub fn imbalance(&self, j: usize) -> Option<BTreeMap<String, f64>> {
        let mut total: BTreeMap<String, f64> = BTreeMap::new();
        for &(i, coef) in &self.reactions[j].stoichiometry {
            let met = &self.metabolites[i];
            let formula = met.formula.as_deref()?;
            for (element, count) in parse_formula(formula)? {
                *total.entry(element).or_insert(0.0) += coef * count;
            }
            if let Some(q) = met.charge {
                *total.entry("charge".to_string()).or_insert(0.0) += coef * f64::from(q);
            }
        }
        total.retain(|_, v| v.abs() > 1e-6);
        Some(total)
    }
}

/// Size and objective summary of a model.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ModelSummary {
    /// Identifier.
    pub id: String,
    /// Name.
    pub name: String,
    /// Reaction count.
    pub reactions: usize,
    /// Metabolite count.
    pub metabolites: usize,
    /// Gene count.
    pub genes: usize,
    /// Compartment count.
    pub compartments: usize,
    /// Exchange reaction count.
    pub exchanges: usize,
    /// Objective reactions and coefficients.
    pub objective: Vec<(String, f64)>,
}

/// Parses a Hill-style formula such as `C6H12O6` into element counts.
/// Formulas with generic groups (`R`, `X`) parse as extra elements; `None` means malformed.
pub fn parse_formula(formula: &str) -> Option<Vec<(String, f64)>> {
    let bytes = formula.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        if !bytes[i].is_ascii_uppercase() {
            return None;
        }
        let start = i;
        i += 1;
        while i < bytes.len() && bytes[i].is_ascii_lowercase() {
            i += 1;
        }
        let element = &formula[start..i];
        let num_start = i;
        while i < bytes.len() && (bytes[i].is_ascii_digit() || bytes[i] == b'.') {
            i += 1;
        }
        let count = if num_start == i { 1.0 } else { formula[num_start..i].parse().ok()? };
        out.push((element.to_string(), count));
    }
    Some(out)
}

/// A value from a kinetic law's parameters, where COBRA files before FBC keep bounds and objectives.
fn local_value(r: &crate::sbml::Reaction, name: &str) -> Option<f64> {
    r.local_parameters.iter().find(|p| p.id == name).and_then(|p| p.value)
}

fn read_reaction(
    r: &sbml::Reaction,
    params: &HashMap<&str, f64>,
    species_index: &HashMap<&str, usize>,
    gene_labels: &HashMap<&str, String>,
    genes: &mut GeneIndex,
    subsystem: &HashMap<&str, &str>,
) -> Result<Reaction> {
    let bound = |name: &Option<String>| -> Result<Option<f64>> {
        match name {
            None => Ok(None),
            Some(p) => match params.get(p.as_str()) {
                Some(v) => Ok(Some(*v)),
                None => Err(Error::Unknown { kind: "bound parameter", id: p.clone() }),
            },
        }
    };
    let lb = bound(&r.lower_bound)?.or_else(|| local_value(r, "LOWER_BOUND")).unwrap_or(if r.reversible {
        f64::NEG_INFINITY
    } else {
        0.0
    });
    let ub = bound(&r.upper_bound)?.or_else(|| local_value(r, "UPPER_BOUND")).unwrap_or(f64::INFINITY);
    let mut merged: BTreeMap<usize, f64> = BTreeMap::new();
    for (refs, sign) in [(&r.reactants, -1.0), (&r.products, 1.0)] {
        for sr in refs {
            let i = *species_index
                .get(sr.species.as_str())
                .ok_or_else(|| Error::Unknown { kind: "species", id: sr.species.clone() })?;
            *merged.entry(i).or_insert(0.0) += sign * sr.stoichiometry.unwrap_or(1.0);
        }
    }
    let gpr = match &r.gene_rule {
        Some(rule) if !rule.trim().is_empty() => {
            let mut intern = |g: &str| match gene_labels.get(g) {
                Some(label) => genes.intern(label),
                None => genes.intern(g.strip_prefix("G_").unwrap_or(g)),
            };
            Gpr::parse(rule, &mut intern)?
        }
        _ => Gpr::None,
    };
    Ok(Reaction {
        id: strip(&r.id, "R_"),
        name: r.name.clone(),
        lb,
        ub,
        stoichiometry: merged.into_iter().filter(|(_, c)| *c != 0.0).collect(),
        gpr,
        subsystem: subsystem.get(r.id.as_str()).map(|s| s.to_string()),
        annotations: r.annotations.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formula_parsing() {
        assert_eq!(parse_formula("C6H12O6").unwrap(), vec![("C".into(), 6.0), ("H".into(), 12.0), ("O".into(), 6.0)]);
        assert_eq!(parse_formula("Fe").unwrap(), vec![("Fe".into(), 1.0)]);
        assert!(parse_formula("c6").is_none());
    }
}
