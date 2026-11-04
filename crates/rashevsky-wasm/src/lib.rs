//! iCH360 flux balance in the browser, on the pure Rust simplex.

use engine::Model;
use engine::lp::simplex::Simplex;
use engine::lp::{Backend, Problem, Sense, Session, Status};
use serde_json::json;
use wasm_bindgen::prelude::*;

const CORE: &[u8] = include_bytes!("../../../models/ich360/iCH360.xml.gz");
const LAYOUT: &str = include_str!("../../../models/ich360/layout.json");

/// The model and a solver session kept between solves, so that each one starts from the last basis.
#[wasm_bindgen]
pub struct Core {
    model: Model,
    session: Box<dyn Session>,
    /// The bounds now in the session, per reaction.
    applied: Vec<(f64, f64)>,
    growth_row: usize,
}

#[wasm_bindgen]
impl Core {
    /// Reads the bundled model.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Result<Core, JsError> {
        let doc = engine::sbml::read::read_bytes(CORE).map_err(|e| JsError::new(&e.to_string()))?;
        let model = Model::from_document(&doc).map_err(|e| JsError::new(&e.to_string()))?;
        let problem = split_problem(&model);
        let growth_row = problem.n_rows() - 1;
        // Far above what this model needs, and short enough that a stuck solve cannot hang the page.
        let solver = Simplex { max_iterations: 50_000, ..Simplex::default() };
        let session = solver.open(&problem).map_err(|e| JsError::new(&e.to_string()))?;
        let applied = model.reactions.iter().map(|r| (r.lb, r.ub)).collect();
        Ok(Core { model, session, applied, growth_row })
    }

    /// The map as JSON.
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
    pub fn solve(&mut self, glucose: f64, oxygen: f64, reactions: &str) -> Result<String, JsError> {
        let fail = |e: engine::Error| JsError::new(&e.to_string());
        let mut wanted: Vec<(f64, f64)> = self.model.reactions.iter().map(|r| (r.lb, r.ub)).collect();
        wanted[self.index("EX_glc__D_e").map_err(fail)?] = (-glucose, 1000.0);
        wanted[self.index("EX_o2_e").map_err(fail)?] = (-oxygen, 1000.0);
        for id in reactions.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            wanted[self.index(id).map_err(fail)?] = (0.0, 0.0);
        }
        match self.least_flux(&wanted).map_err(fail)? {
            Some((growth, fluxes)) => Ok(json!({"status": "optimal", "growth": growth, "fluxes": fluxes}).to_string()),
            None => Ok(json!({"status": "none", "growth": null, "fluxes": vec![0.0; wanted.len()]}).to_string()),
        }
    }
}

/// Reactions split into forward and reverse columns, whose sum is the total flux, and a last row that
/// will hold growth at its optimum.
fn split_problem(m: &Model) -> Problem {
    let n = m.reactions.len();
    let mut p = Problem { sense: Some(Sense::Maximize), ..Problem::default() };
    for r in &m.reactions {
        p.add_col(0.0, r.lb.max(0.0), r.ub.max(0.0), false);
    }
    for r in &m.reactions {
        p.add_col(0.0, (-r.ub).max(0.0), (-r.lb).max(0.0), false);
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
    p.add_row(f64::NEG_INFINITY, f64::INFINITY, objective);
    p
}

impl Core {
    fn index(&self, id: &str) -> engine::Result<usize> {
        self.model.reaction_index(id)
    }

    /// The two steps of pFBA: the most growth, then the least total flux that keeps it. `None` when the
    /// bounds leave no steady state, and an error when the solver gives no verdict.
    fn least_flux(&mut self, wanted: &[(f64, f64)]) -> engine::Result<Option<(f64, Vec<f64>)>> {
        let n = wanted.len();
        for (j, &(lb, ub)) in wanted.iter().enumerate() {
            if self.applied[j] != (lb, ub) {
                self.session.set_col_bounds(j, lb.max(0.0), ub.max(0.0))?;
                self.session.set_col_bounds(n + j, (-ub).max(0.0), (-lb).max(0.0))?;
                self.applied[j] = (lb, ub);
            }
        }
        let mut growth = vec![0.0; 2 * n];
        for &(j, c) in &self.model.objective {
            growth[j] = c;
            growth[n + j] = -c;
        }
        self.session.set_row_bounds(self.growth_row, f64::NEG_INFINITY, f64::INFINITY)?;
        self.session.set_objective(Sense::Maximize, &growth)?;
        let first = self.session.solve()?;
        match first.status {
            Status::Optimal => {}
            Status::Infeasible => return Ok(None),
            other => return Err(engine::Error::Solver(format!("growth step ended with {other:?}"))),
        }
        self.session.set_row_bounds(self.growth_row, first.objective, f64::INFINITY)?;
        self.session.set_objective(Sense::Minimize, &vec![1.0; 2 * n])?;
        let s = self.session.solve()?;
        if s.status != Status::Optimal {
            return Err(engine::Error::Solver(format!("least flux step ended with {:?}", s.status)));
        }
        Ok(Some((first.objective, (0..n).map(|j| s.x[j] - s.x[n + j]).collect())))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use engine::cb;

    fn value(text: &str) -> serde_json::Value {
        serde_json::from_str(text).unwrap()
    }

    fn core() -> Core {
        Core::new().map_err(|_| "new").unwrap()
    }

    #[test]
    fn matches_the_reference_growth() {
        let out = value(&core().solve(10.0, 1000.0, "").map_err(|_| "solve").unwrap());
        assert!((out["growth"].as_f64().unwrap() - 0.8656270306509).abs() < 1e-7);
    }

    #[test]
    fn a_knockout_changes_growth_as_in_the_native_engine() {
        let mut core = core();
        let out = value(&core.solve(10.0, 1000.0, "PGI").map_err(|_| "solve").unwrap());
        let mut m = core.model.clone();
        m.set_bounds("PGI", 0.0, 0.0).unwrap();
        let native = cb::fba_with(&m, &Simplex::default()).unwrap();
        assert!((out["growth"].as_f64().unwrap() - native.objective).abs() < 1e-8);
    }

    #[test]
    fn anaerobic_growth_is_slower() {
        let out = value(&core().solve(10.0, 0.0, "").map_err(|_| "solve").unwrap());
        assert!((out["growth"].as_f64().unwrap() - 0.1565800).abs() < 1e-6);
    }

    /// A long session must give what a new one gives at every step, since the page keeps one open.
    #[test]
    fn a_long_session_matches_a_fresh_solve_at_every_step() {
        let mut long = core();
        let knockouts = ["", "PGI", "PGI,PFK", "ENO", "CS", "", "ATPS4rpp", "GLNS,PYK", "CYTBO3_4pp", ""];
        let mut no_state = 0;
        for step in 0..60 {
            let glucose = [10.0, 2.5, 0.0, 7.0, 20.0][step % 5];
            let oxygen = [1000.0, 0.0, 12.0, 3.0][step % 4];
            let cut = knockouts[step % knockouts.len()];
            let kept = value(&long.solve(glucose, oxygen, cut).map_err(|_| "solve").unwrap());
            let fresh = value(&core().solve(glucose, oxygen, cut).map_err(|_| "solve").unwrap());
            assert_eq!(kept["status"], fresh["status"], "step {step}");
            if kept["status"] == "none" {
                no_state += 1;
                continue;
            }
            let (a, b) = (kept["growth"].as_f64().unwrap(), fresh["growth"].as_f64().unwrap());
            assert!((a - b).abs() < 1e-7, "step {step}: {a} {b}");
            let total = |o: &serde_json::Value| -> f64 {
                o["fluxes"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap().abs()).sum()
            };
            assert!((total(&kept) - total(&fresh)).abs() < 1e-5, "step {step}: {} {}", total(&kept), total(&fresh));
            // The fluxes themselves are a steady state inside the bounds asked for.
            let fluxes: Vec<f64> = kept["fluxes"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
            for row in long.model.to_lp().rows {
                let balance: f64 = row.iter().map(|&(j, c)| c * fluxes[j]).sum();
                assert!(balance.abs() < 1e-6, "step {step}: {balance}");
            }
            for id in cut.split(',').filter(|s| !s.is_empty()) {
                assert_eq!(fluxes[long.model.reaction_index(id).unwrap()], 0.0, "step {step}, {id}");
            }
        }
        assert!(no_state > 0 && no_state < 60, "{no_state} of the steps had no steady state");
    }

    /// Compares a long session with a new one wherever the session reports no steady state, and on every
    /// fiftieth step otherwise, since a new session costs far more than a kept one.
    fn walk(moves: &[(f64, f64, String)]) {
        let mut long = core();
        for (step, (glucose, oxygen, cut)) in moves.iter().enumerate() {
            let kept =
                long.solve(*glucose, *oxygen, cut).map_err(|_| format!("step {step} threw: {glucose} {oxygen} {cut}"));
            let kept = value(&kept.unwrap());
            if kept["status"] != "optimal" || step % 50 == 0 {
                let fresh = value(&core().solve(*glucose, *oxygen, cut).map_err(|_| "fresh").unwrap());
                assert_eq!(
                    kept["status"], fresh["status"],
                    "step {step}: {glucose} {oxygen} {cut}: {}",
                    fresh["growth"]
                );
                if let (Some(a), Some(b)) = (kept["growth"].as_f64(), fresh["growth"].as_f64()) {
                    assert!((a - b).abs() < 1e-6, "step {step}: {a} {b}");
                }
            }
        }
    }

    /// A small generator, so that the walks repeat exactly.
    fn mulberry(state: &mut u32) -> f64 {
        *state = state.wrapping_add(0x6D2B79F5);
        let mut t = *state;
        t = (t ^ (t >> 15)).wrapping_mul(t | 1);
        t ^= t.wrapping_add((t ^ (t >> 7)).wrapping_mul(t | 61));
        f64::from(t ^ (t >> 14)) / 4294967296.0
    }

    /// Dragging the sliders for a long time. Ended in a false verdict at step 3767 when a kept basis was
    /// believed without a check.
    #[test]
    fn long_slider_walks_keep_the_right_verdict() {
        let mut state = 21u32;
        let (mut glucose, mut oxygen) = (10.0, 60.0);
        let mut moves = Vec::new();
        for _ in 0..4000 {
            if mulberry(&mut state) < 0.6 {
                glucose = (mulberry(&mut state) * 61.0).floor() * 0.5;
            } else {
                oxygen = (mulberry(&mut state) * 61.0).floor();
            }
            moves.push((glucose, if oxygen >= 60.0 { 1000.0 } else { oxygen }, String::new()));
        }
        walk(&moves);
    }

    /// Deleting each reaction in turn and restoring it. Ended in a false verdict at step 3327 when a kept
    /// basis was believed without a check.
    #[test]
    fn deleting_every_reaction_in_turn_keeps_the_right_verdict() {
        let ids: Vec<String> =
            core().model.reactions.iter().map(|r| r.id.clone()).filter(|id| id != "Biomass").collect();
        let mut moves = Vec::new();
        for round in 0..8 {
            let (glucose, oxygen) = ([10.0, 7.5, 12.0, 20.0][round % 4], [1000.0, 20.0, 1000.0, 8.0][round % 4]);
            for id in &ids {
                moves.push((glucose, oxygen, id.clone()));
                moves.push((glucose, oxygen, String::new()));
            }
        }
        walk(&moves);
    }

    #[test]
    fn an_unknown_reaction_changes_nothing_in_the_session() {
        let mut core = core();
        let before = value(&core.solve(10.0, 1000.0, "").map_err(|_| "solve").unwrap());
        assert!(core.index("NOPE").is_err());
        let after = value(&core.solve(10.0, 1000.0, "").map_err(|_| "solve").unwrap());
        assert!((before["growth"].as_f64().unwrap() - after["growth"].as_f64().unwrap()).abs() < 1e-9);
    }
}
