//! `e_coli_core` flux balance in the browser, on the pure Rust simplex.

use engine::Model;
use engine::cb;
use engine::lp::simplex::Simplex;
use engine::lp::{Backend, Problem, Sense, Status};
use serde_json::json;
use wasm_bindgen::prelude::*;

const CORE: &[u8] = include_bytes!("../../../models/e_coli_core/e_coli_core.xml.gz");
const LAYOUT: &str = include_str!("../../../models/e_coli_core/layout.json");

/// The model and the settings of the demo.
#[wasm_bindgen]
pub struct Core {
    model: Model,
    backend: Simplex,
}

#[wasm_bindgen]
impl Core {
    /// Reads the bundled model.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Result<Core, JsError> {
        let doc = engine::sbml::read::read_bytes(CORE).map_err(|e| JsError::new(&e.to_string()))?;
        let model = Model::from_document(&doc).map_err(|e| JsError::new(&e.to_string()))?;
        Ok(Core { model, backend: Simplex::default() })
    }

    /// The hand-placed map as JSON.
    pub fn layout(&self) -> String {
        LAYOUT.to_string()
    }

    /// Reactions and their genes as JSON.
    pub fn reactions(&self) -> String {
        let list: Vec<_> = self
            .model
            .reactions
            .iter()
            .map(|r| {
                let mut genes = Vec::new();
                r.gpr.genes(&mut genes);
                json!({
                    "id": r.id,
                    "name": r.name,
                    "rule": r.gpr.to_rule(&|g| self.model.genes[g].id.clone()),
                    "genes": genes.iter().map(|&g| self.model.genes[g].id.clone()).collect::<Vec<_>>(),
                })
            })
            .collect();
        json!(list).to_string()
    }

    /// Parsimonious FBA with the given uptake limits and deletions; `reactions` holds identifiers
    /// separated by commas. Returns `{status, growth, fluxes}` as JSON.
    pub fn solve(&self, glucose: f64, oxygen: f64, reactions: &str) -> Result<String, JsError> {
        let fail = |e: engine::Error| JsError::new(&e.to_string());
        let mut m = self.model.clone();
        m.set_bounds("EX_glc__D_e", -glucose, 1000.0).map_err(fail)?;
        m.set_bounds("EX_o2_e", -oxygen, 1000.0).map_err(fail)?;
        for id in reactions.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            m.set_bounds(id, 0.0, 0.0).map_err(fail)?;
        }
        let first = cb::fba_with(&m, &self.backend).map_err(fail)?;
        if first.status != Status::Optimal {
            return Ok(json!({"status": "none", "growth": null, "fluxes": vec![0.0; m.reactions.len()]}).to_string());
        }
        let fluxes = self.least_flux(&m, first.objective).map_err(fail)?;
        Ok(json!({"status": "optimal", "growth": first.objective, "fluxes": fluxes}).to_string())
    }
}

impl Core {
    /// The second step of pFBA, with reactions split into forward and reverse parts.
    fn least_flux(&self, m: &Model, growth: f64) -> engine::Result<Vec<f64>> {
        let n = m.reactions.len();
        let mut p = Problem { sense: Some(Sense::Minimize), ..Problem::default() };
        for r in &m.reactions {
            p.add_col(1.0, r.lb.max(0.0), r.ub.max(0.0), false);
        }
        for r in &m.reactions {
            p.add_col(1.0, (-r.ub).max(0.0), (-r.lb).max(0.0), false);
        }
        for row in m.to_lp().rows {
            let mut coefs = Vec::with_capacity(row.len() * 2);
            for (j, c) in row {
                coefs.push((j, c));
                coefs.push((n + j, -c));
            }
            p.add_row(0.0, 0.0, coefs);
        }
        let mut objective = Vec::new();
        for &(j, c) in &m.objective {
            objective.push((j, c));
            objective.push((n + j, -c));
        }
        p.add_row(growth, f64::INFINITY, objective);
        let mut session = self.backend.open(&p)?;
        let s = session.solve()?;
        if s.status != Status::Optimal {
            return Err(engine::Error::Solver(format!("least flux step ended with {:?}", s.status)));
        }
        Ok((0..n).map(|j| s.x[j] - s.x[n + j]).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn value(text: &str) -> serde_json::Value {
        serde_json::from_str(text).unwrap()
    }

    #[test]
    fn matches_the_reference_growth() {
        let core = Core::new().map_err(|_| "new").unwrap();
        let out = value(&core.solve(10.0, 1000.0, "").map_err(|_| "solve").unwrap());
        assert!((out["growth"].as_f64().unwrap() - 0.8739215069684).abs() < 1e-7);
    }

    #[test]
    fn a_knockout_changes_growth_as_in_the_native_engine() {
        let core = Core::new().map_err(|_| "new").unwrap();
        let out = value(&core.solve(10.0, 1000.0, "PGI").map_err(|_| "solve").unwrap());
        let mut m = core.model.clone();
        m.set_bounds("PGI", 0.0, 0.0).unwrap();
        let native = cb::fba_with(&m, &core.backend).unwrap();
        assert!((out["growth"].as_f64().unwrap() - native.objective).abs() < 1e-8);
    }

    #[test]
    fn anaerobic_growth_is_slower() {
        let core = Core::new().map_err(|_| "new").unwrap();
        let out = value(&core.solve(10.0, 0.0, "").map_err(|_| "solve").unwrap());
        assert!((out["growth"].as_f64().unwrap() - 0.21166).abs() < 1e-4);
    }
}
