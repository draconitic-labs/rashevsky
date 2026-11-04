//! Two-dimensional layout of a metabolic network, by a seeded force-directed method.
//!
//! The graph joins reactions to their primary metabolites. Currency metabolites (those joined to many
//! reactions, such as water and ATP) stay out of the layout, or a map would collapse around them.
//! This is the Fruchterman-Reingold (1991) method with a cooling schedule; it draws a readable network
//! and has no biological meaning, so curated layouts are preferred where they exist.

use crate::model::Model;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};

/// A placed metabolite or reaction.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Placed {
    /// Index into the model's metabolites or reactions.
    pub index: usize,
    /// Horizontal position.
    pub x: f64,
    /// Vertical position.
    pub y: f64,
}

/// A layout of part of a model.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Layout {
    /// Primary metabolites.
    pub metabolites: Vec<Placed>,
    /// Reactions, at their midpoints.
    pub reactions: Vec<Placed>,
    /// Width of the drawing.
    pub width: f64,
    /// Height of the drawing.
    pub height: f64,
    /// Metabolites treated as currency in this layout.
    pub currency: Vec<usize>,
}

/// Metabolites joined to at least `threshold` of the reactions are treated as currency.
pub fn currency_metabolites(model: &Model, reactions: &[usize], threshold: usize) -> Vec<usize> {
    let mut degree = vec![0usize; model.metabolites.len()];
    for &j in reactions {
        for &(i, _) in &model.reactions[j].stoichiometry {
            degree[i] += 1;
        }
    }
    (0..degree.len()).filter(|&i| degree[i] >= threshold).collect()
}

/// Lays out the given reactions.
pub fn force_layout(model: &Model, reactions: &[usize], seed: u64, iterations: usize) -> Layout {
    let threshold = (reactions.len() / 8).max(6);
    let currency = currency_metabolites(model, reactions, threshold);
    let is_currency = |i: usize| currency.binary_search(&i).is_ok();
    let mut met_node = vec![usize::MAX; model.metabolites.len()];
    let mut mets = Vec::new();
    for &j in reactions {
        for &(i, _) in &model.reactions[j].stoichiometry {
            if !is_currency(i) && met_node[i] == usize::MAX {
                met_node[i] = mets.len();
                mets.push(i);
            }
        }
    }
    let n_m = mets.len();
    let n = n_m + reactions.len();
    let mut edges = Vec::new();
    for (k, &j) in reactions.iter().enumerate() {
        for &(i, _) in &model.reactions[j].stoichiometry {
            if !is_currency(i) {
                edges.push((met_node[i], n_m + k));
            }
        }
    }
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let side = (n as f64).sqrt().max(1.0) * 120.0;
    let mut pos: Vec<(f64, f64)> = (0..n).map(|_| (rng.random::<f64>() * side, rng.random::<f64>() * side)).collect();
    let k = side / (n as f64).sqrt().max(1.0);
    let mut temperature = side / 6.0;
    let cooling = (0.02f64).powf(1.0 / iterations.max(1) as f64);
    for _ in 0..iterations {
        let mut disp = vec![(0.0f64, 0.0f64); n];
        for a in 0..n {
            for b in (a + 1)..n {
                let (dx, dy) = (pos[a].0 - pos[b].0, pos[a].1 - pos[b].1);
                let d2 = (dx * dx + dy * dy).max(1e-3);
                let force = k * k / d2;
                disp[a].0 += dx * force;
                disp[a].1 += dy * force;
                disp[b].0 -= dx * force;
                disp[b].1 -= dy * force;
            }
        }
        for &(a, b) in &edges {
            let (dx, dy) = (pos[a].0 - pos[b].0, pos[a].1 - pos[b].1);
            let d = (dx * dx + dy * dy).sqrt().max(1e-3);
            let force = d / k;
            disp[a].0 -= dx * force;
            disp[a].1 -= dy * force;
            disp[b].0 += dx * force;
            disp[b].1 += dy * force;
        }
        // Gravity toward the centre keeps separate components of the network from drifting apart.
        for a in 0..n {
            disp[a].0 -= (pos[a].0 - side / 2.0) * 0.35;
            disp[a].1 -= (pos[a].1 - side / 2.0) * 0.35;
        }
        for a in 0..n {
            let (dx, dy) = disp[a];
            let d = (dx * dx + dy * dy).sqrt().max(1e-3);
            let step = d.min(temperature);
            pos[a].0 += dx / d * step;
            pos[a].1 += dy / d * step;
        }
        temperature *= cooling;
    }
    let min_x = pos.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
    let min_y = pos.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
    let max_x = pos.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
    let max_y = pos.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
    let margin = 60.0;
    let place = |p: (f64, f64)| (p.0 - min_x + margin, p.1 - min_y + margin);
    Layout {
        metabolites: mets
            .iter()
            .enumerate()
            .map(|(a, &i)| {
                let (x, y) = place(pos[a]);
                Placed { index: i, x, y }
            })
            .collect(),
        reactions: reactions
            .iter()
            .enumerate()
            .map(|(a, &j)| {
                let (x, y) = place(pos[n_m + a]);
                Placed { index: j, x, y }
            })
            .collect(),
        width: (max_x - min_x).max(1.0) + 2.0 * margin,
        height: (max_y - min_y).max(1.0) + 2.0 * margin,
        currency,
    }
}
