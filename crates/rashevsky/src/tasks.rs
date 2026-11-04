//! Metabolic task checks in the task-list format of RAVEN and Human-GEM.
//!
//! A task states allowed inputs and outputs, extra reactions and bound changes; it passes when a
//! steady-state flux satisfying them exists (or does not, for tasks marked `SHOULD FAIL`).
//! Inputs and outputs act on the balance row of a metabolite: an input allows net consumption
//! `S v in [-UB, -LB]`, an output allows net production `S v <= UB` and requires `S v >= LB` when
//! `LB > 0`. Semantics follow `checkTasks` of RAVEN and its Python port, with boundary reactions
//! closed so the task's inputs and outputs are the complete boundary, as in Human-GEM's own test.

use crate::error::{Error, Result};
use crate::lp::{self, Sense, Status};
use crate::model::Model;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// One metabolic task. Bounds are `(name, lower, upper)`.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Task {
    /// Task identifier.
    pub id: String,
    /// Description.
    pub description: String,
    /// True when the task must be infeasible.
    pub should_fail: bool,
    /// Inputs as `name[compartment]`, or `ALLMETS` and `ALLMETSIN[compartment]`.
    pub inputs: Vec<(String, f64, f64)>,
    /// Outputs, in the same form.
    pub outputs: Vec<(String, f64, f64)>,
    /// Extra reactions as equations.
    pub equations: Vec<(String, f64, f64)>,
    /// Reactions whose bounds change.
    pub changed: Vec<(String, f64, f64)>,
}

/// Outcome of one task.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TaskResult {
    /// Task identifier.
    pub id: String,
    /// Description.
    pub description: String,
    /// Verdict, accounting for `should_fail`.
    pub passed: bool,
    /// Whether a steady-state flux exists.
    pub feasible: bool,
    /// Reason the task could not be evaluated.
    pub error: Option<String>,
}

fn number(text: &str, default: f64) -> Result<f64> {
    let t = text.trim();
    if t.is_empty() {
        return Ok(default);
    }
    t.parse().map_err(|_| Error::Invalid(format!("bad number `{t}` in task list")))
}

/// Parses a tab-separated task list.
pub fn parse_task_list(text: &str) -> Result<Vec<Task>> {
    let rows: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    let header_at = rows
        .iter()
        .position(|r| r.iter().any(|c| c.trim().eq_ignore_ascii_case("ID")))
        .ok_or_else(|| Error::Invalid("task list has no ID column".into()))?;
    let header: Vec<String> = rows[header_at].iter().map(|c| c.trim().to_uppercase()).collect();
    let col: HashMap<&str, usize> = header.iter().enumerate().map(|(i, c)| (c.as_str(), i)).collect();
    let id_col = col["ID"];
    let cell = |row: &[&str], name: &str| -> String {
        col.get(name).and_then(|&i| row.get(i)).map_or(String::new(), |c| c.trim().to_string())
    };
    let split = |text: &str| -> Vec<String> {
        text.split(';').map(str::trim).filter(|s| !s.is_empty()).map(str::to_string).collect()
    };
    let mut tasks: Vec<Task> = Vec::new();
    for row in &rows[header_at + 1..] {
        if row.iter().all(|c| c.trim().is_empty()) {
            continue;
        }
        // A non-empty cell before the ID column marks a comment row.
        if row.iter().take(id_col).any(|c| !c.trim().is_empty()) {
            continue;
        }
        let id = cell(row, "ID");
        if id.starts_with('#') {
            continue;
        }
        if !id.is_empty() {
            let fail = cell(row, "SHOULD FAIL").to_lowercase();
            tasks.push(Task {
                id,
                description: cell(row, "DESCRIPTION"),
                should_fail: !matches!(fail.as_str(), "" | "0" | "false" | "no"),
                ..Task::default()
            });
        }
        let Some(task) = tasks.last_mut() else { continue };
        let input = cell(row, "IN");
        if !input.is_empty() {
            let (lb, ub) = (number(&cell(row, "IN LB"), 0.0)?, number(&cell(row, "IN UB"), 1000.0)?);
            task.inputs.extend(split(&input).into_iter().map(|m| (m, lb, ub)));
        }
        let output = cell(row, "OUT");
        if !output.is_empty() {
            let (lb, ub) = (number(&cell(row, "OUT LB"), 0.0)?, number(&cell(row, "OUT UB"), 1000.0)?);
            task.outputs.extend(split(&output).into_iter().map(|m| (m, lb, ub)));
        }
        let equation = cell(row, "EQU");
        if !equation.is_empty() {
            let default_lb = if equation.contains("<=>") { -1000.0 } else { 0.0 };
            let lb = number(&cell(row, "EQU LB"), default_lb)?;
            let ub = number(&cell(row, "EQU UB"), 1000.0)?;
            task.equations.push((equation, lb, ub));
        }
        let changed = cell(row, "CHANGED RXN");
        if !changed.is_empty() {
            let lb = number(&cell(row, "CHANGED LB"), -1000.0)?;
            let ub = number(&cell(row, "CHANGED UB"), 1000.0)?;
            task.changed.extend(split(&changed).into_iter().map(|r| (r, lb, ub)));
        }
    }
    Ok(tasks)
}

/// Lookup of metabolites by `NAME[COMP]` (several may share one) and by compartment.
struct Names {
    by_name: HashMap<String, Vec<usize>>,
    by_compartment: HashMap<String, Vec<usize>>,
}

impl Names {
    fn new(model: &Model) -> Self {
        let mut by_name: HashMap<String, Vec<usize>> = HashMap::new();
        let mut by_compartment: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, m) in model.metabolites.iter().enumerate() {
            by_name.entry(format!("{}[{}]", m.name, m.compartment).to_uppercase()).or_default().push(i);
            by_compartment.entry(m.compartment.to_uppercase()).or_default().push(i);
        }
        Self { by_name, by_compartment }
    }
}

/// Splits `2 ATP[c] + H2O[c]` into `(coefficient, name)`.
fn side(text: &str) -> Vec<(f64, String)> {
    text.split(" + ")
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|term| match term.split_once(' ') {
            Some((head, rest)) => match head.parse::<f64>() {
                Ok(c) => (c, rest.trim().to_string()),
                Err(_) => (1.0, term.to_string()),
            },
            None => (1.0, term.to_string()),
        })
        .collect()
}

fn parse_equation(text: &str) -> Result<Vec<(f64, String)>> {
    for arrow in ["<=>", "-->", "->", "=>"] {
        if let Some((left, right)) = text.split_once(arrow) {
            let mut out: Vec<(f64, String)> = side(left).into_iter().map(|(c, m)| (-c, m)).collect();
            out.extend(side(right));
            return Ok(out);
        }
    }
    Err(Error::Invalid(format!("no reaction arrow in `{text}`")))
}

fn evaluate(
    model: &Model,
    base: &lp::Problem,
    names: &Names,
    task: &Task,
    backend: &dyn lp::Backend,
) -> Result<std::result::Result<bool, String>> {
    let mut p = base.clone();
    let mut bounds: HashMap<usize, (f64, f64)> = HashMap::new();
    let mut forced_in = std::collections::HashSet::new();
    let mut forced_out = std::collections::HashSet::new();
    for (entries, is_input) in [(&task.inputs, true), (&task.outputs, false)] {
        let classify = |token: &str| -> (u8, String) {
            let upper = token.to_uppercase();
            if upper == "ALLMETS" {
                (0, String::new())
            } else if let Some(rest) = upper.strip_prefix("ALLMETSIN[").and_then(|r| r.strip_suffix(']')) {
                (1, rest.to_string())
            } else {
                (2, upper)
            }
        };
        // Bulk tokens first so that specific metabolites override them.
        let mut ordered: Vec<&(String, f64, f64)> = entries.iter().collect();
        ordered.sort_by_key(|(t, _, _)| classify(t).0 == 2);
        for (token, lb, ub) in ordered {
            let (kind, arg) = classify(token);
            let ids: Vec<usize> = match kind {
                0 => (0..model.metabolites.len()).collect(),
                1 => names.by_compartment.get(&arg).cloned().unwrap_or_default(),
                _ => match names.by_name.get(&arg) {
                    Some(v) => v.clone(),
                    None => return Ok(Err(format!("unknown metabolite {token}"))),
                },
            };
            for i in ids {
                // As RAVEN: one metabolite cannot be both forced in and forced out.
                if kind == 2 && *lb > 0.0 {
                    let (this, other) =
                        if is_input { (&mut forced_in, &forced_out) } else { (&mut forced_out, &forced_in) };
                    if other.contains(&i) {
                        return Ok(Err(format!("{token} has a nonzero lower bound both as input and as output")));
                    }
                    this.insert(i);
                }
                let b = bounds.entry(i).or_insert((0.0, 0.0));
                if is_input {
                    b.0 = -ub;
                    if kind == 2 {
                        b.1 = -lb;
                    }
                } else {
                    b.1 = *ub;
                    if kind == 2 && *lb > 0.0 {
                        b.0 = *lb;
                    }
                }
            }
        }
    }
    for (&i, &(lb, ub)) in &bounds {
        if lb > ub {
            return Ok(Err(format!("the limits on {} contradict each other", model.metabolites[i].id)));
        }
        p.row_lb[i] = lb;
        p.row_ub[i] = ub;
    }
    // Extra reactions. Metabolites not in the model become new balanced rows.
    let mut extra: HashMap<String, usize> = HashMap::new();
    for (equation, lb, ub) in &task.equations {
        let terms = parse_equation(equation)?;
        let col = p.add_col(0.0, *lb, *ub, false);
        for (coef, name) in terms {
            let key = name.to_uppercase();
            let rows: Vec<usize> = match names.by_name.get(&key) {
                Some(v) => vec![v[0]],
                None => {
                    let row = *extra.entry(key).or_insert_with(|| p.add_row(0.0, 0.0, Vec::new()));
                    vec![row]
                }
            };
            // A metabolite on both sides of an equation is one coefficient.
            for row in rows {
                match p.rows[row].iter_mut().find(|(c, _)| *c == col) {
                    Some(entry) => entry.1 += coef,
                    None => p.rows[row].push((col, coef)),
                }
            }
        }
    }
    for (rid, lb, ub) in &task.changed {
        let j = match model.reaction_index(rid) {
            Ok(j) => j,
            Err(e) => return Ok(Err(e.to_string())),
        };
        p.col_lb[j] = *lb;
        p.col_ub[j] = *ub;
    }
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    Ok(match s.status {
        Status::Optimal => Ok(true),
        Status::Infeasible => Ok(false),
        other => Err(format!("solver ended with {other:?}")),
    })
}

/// Evaluates a task list. With `close_boundaries`, reactions with one metabolite of unit coefficient
/// are closed, so the task's inputs and outputs define the whole boundary of the system.
pub fn check_tasks(model: &Model, tasks: &[Task], close_boundaries: bool) -> Result<Vec<TaskResult>> {
    let mut base = model.to_lp();
    base.obj.iter_mut().for_each(|c| *c = 0.0);
    base.sense = Some(Sense::Maximize);
    if close_boundaries {
        for j in 0..model.reactions.len() {
            if model.exchange_coefficient(j).is_some_and(|c| c.abs() == 1.0) {
                base.col_lb[j] = 0.0;
                base.col_ub[j] = 0.0;
            }
        }
    }
    let names = Names::new(model);
    let backend = lp::default_backend();
    let parts: Vec<Result<TaskResult>> = crate::parallel::map_chunks(tasks, 1, |slice| {
        let task = &slice[0];
        let outcome = evaluate(model, &base, &names, task, backend.as_ref())?;
        Ok(match outcome {
            Ok(feasible) => TaskResult {
                id: task.id.clone(),
                description: task.description.clone(),
                passed: feasible != task.should_fail,
                feasible,
                error: None,
            },
            Err(reason) => TaskResult {
                id: task.id.clone(),
                description: task.description.clone(),
                passed: false,
                feasible: false,
                error: Some(reason),
            },
        })
    });
    parts.into_iter().collect()
}
