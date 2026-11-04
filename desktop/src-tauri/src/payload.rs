//! What the interface receives about a model.

use rashevsky::Model;
use rashevsky::io::layout::force_layout;
use serde::Serialize;
use serde_json::{Map, Value, json};
use std::collections::BTreeSet;

/// A reaction as the interface shows it.
#[derive(Serialize)]
pub struct ReactionView {
    pub id: String,
    pub name: String,
    pub lb: f64,
    pub ub: f64,
    pub equation: String,
    pub rule: String,
    pub genes: Vec<String>,
    pub subsystem: Option<String>,
    pub boundary: bool,
}

/// A model as the interface shows it.
#[derive(Serialize)]
pub struct ModelView {
    pub id: String,
    pub name: String,
    pub reactions: Vec<ReactionView>,
    pub metabolites: Vec<Value>,
    pub genes: Vec<Value>,
    pub objective: Vec<String>,
    pub compartments: usize,
    pub subsystems: Vec<String>,
    pub medium: Vec<Value>,
    pub curated_layout: bool,
    pub balance: Value,
}

fn side(model: &Model, stoich: &[(usize, f64)], sign: f64) -> String {
    let mut parts = Vec::new();
    for &(i, c) in stoich {
        if c * sign > 0.0 {
            let name = &model.metabolites[i].id;
            let amount = c.abs();
            if (amount - 1.0).abs() < 1e-9 {
                parts.push(name.clone());
            } else {
                parts.push(format!("{amount:.4} {name}").replace(".0000", ""));
            }
        }
    }
    parts.join(" + ")
}

pub fn model_view(model: &Model, curated_layout: bool) -> ModelView {
    let reactions = model
        .reactions
        .iter()
        .map(|r| {
            let mut gene_ids = Vec::new();
            r.gpr.genes(&mut gene_ids);
            gene_ids.sort_unstable();
            gene_ids.dedup();
            let arrow = if r.lb < 0.0 && r.ub > 0.0 {
                "⇌"
            } else if r.ub <= 0.0 && r.lb < 0.0 {
                "←"
            } else {
                "→"
            };
            ReactionView {
                id: r.id.clone(),
                name: r.name.clone(),
                lb: r.lb,
                ub: r.ub,
                equation: format!("{} {arrow} {}", side(model, &r.stoichiometry, -1.0), side(model, &r.stoichiometry, 1.0)),
                rule: r.gpr.to_rule(&|g| model.genes[g].id.clone()),
                genes: gene_ids.iter().map(|&g| model.genes[g].id.clone()).collect(),
                subsystem: r.subsystem.clone(),
                boundary: r.is_boundary(),
            }
        })
        .collect();
    let subsystems: BTreeSet<String> = model.reactions.iter().filter_map(|r| r.subsystem.clone()).collect();
    let mut medium: Vec<Value> = rashevsky::cb::medium(model)
        .into_iter()
        .filter(|(_, limit)| limit.is_finite())
        .map(|(id, limit)| json!({"id": id, "limit": limit}))
        .collect();
    medium.sort_by(|a, b| b["limit"].as_f64().partial_cmp(&a["limit"].as_f64()).unwrap_or(std::cmp::Ordering::Equal));
    let (mut balanced, mut unknown, mut imbalanced) = (0, 0, 0);
    for (j, r) in model.reactions.iter().enumerate() {
        if r.is_boundary() {
            continue;
        }
        match model.imbalance(j) {
            None => unknown += 1,
            Some(map) if map.is_empty() => balanced += 1,
            Some(_) => imbalanced += 1,
        }
    }
    ModelView {
        id: model.id.clone(),
        name: model.name.clone(),
        reactions,
        metabolites: model
            .metabolites
            .iter()
            .map(|m| json!({"id": m.id, "name": m.name, "compartment": m.compartment, "formula": m.formula, "charge": m.charge}))
            .collect(),
        genes: model.genes.iter().map(|g| json!({"id": g.id, "name": g.name})).collect(),
        objective: model.objective.iter().map(|&(j, _)| model.reactions[j].id.clone()).collect(),
        compartments: model.compartments.len(),
        subsystems: subsystems.into_iter().collect(),
        medium,
        curated_layout,
        balance: json!({"balanced": balanced, "formula_missing": unknown, "imbalanced": imbalanced}),
    }
}

/// A force-directed map of the given reactions, in the form of the hand-placed layout files.
pub fn computed_layout(model: &Model, reactions: &[usize]) -> Value {
    let layout = force_layout(model, reactions, 1, if reactions.len() > 200 { 220 } else { 400 });
    let mut nodes = Map::new();
    for p in &layout.metabolites {
        let m = &model.metabolites[p.index];
        let mut label = if m.name.is_empty() { m.id.clone() } else { m.name.clone() };
        if label.chars().count() > 18 {
            label = label.chars().take(17).collect::<String>() + "…";
        }
        nodes.insert(m.id.clone(), json!({"x": p.x, "y": p.y, "label": label, "external": m.compartment == "e"}));
    }
    let drawn: std::collections::HashSet<usize> = layout.metabolites.iter().map(|p| p.index).collect();
    let mut entries = Map::new();
    for p in &layout.reactions {
        let r = &model.reactions[p.index];
        let list = |sign: f64| -> Vec<&str> {
            r.stoichiometry
                .iter()
                .filter(|(i, c)| c * sign > 0.0 && drawn.contains(i))
                .map(|&(i, _)| model.metabolites[i].id.as_str())
                .collect()
        };
        entries.insert(r.id.clone(), json!({"s": list(-1.0), "p": list(1.0), "m": [p.x, p.y], "boundary": r.is_boundary()}));
    }
    json!({"world": [layout.width, layout.height], "nodes": nodes, "regions": [], "reactions": entries, "source": "force-directed"})
}
