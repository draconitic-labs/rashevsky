//! The growth medium: exchange uptake limits, and minimal media that sustain a required growth.

use crate::error::{Error, Result};
use crate::lp::{self, Problem, Sense};
use crate::model::Model;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Uptake limit per exchange reaction that can take up metabolite, as a positive number.
///
/// Exchange reactions follow the BiGG convention `met <->`, so uptake is negative flux.
pub fn medium(model: &Model) -> BTreeMap<String, f64> {
    let mut out = BTreeMap::new();
    for j in model.exchanges() {
        let r = &model.reactions[j];
        let coef = model.exchange_coefficient(j).unwrap_or(-1.0);
        let limit = if coef < 0.0 { -r.lb } else { r.ub };
        if limit > 0.0 {
            out.insert(r.id.clone(), limit);
        }
    }
    out
}

/// Sets the medium: listed exchanges get the given uptake limit, all others are closed for uptake.
pub fn set_medium(model: &mut Model, uptake: &BTreeMap<String, f64>) -> Result<()> {
    for id in uptake.keys() {
        model.reaction_index(id)?;
    }
    for j in model.exchanges() {
        let coef = model.exchange_coefficient(j).unwrap_or(-1.0);
        let r = &mut model.reactions[j];
        let limit = uptake.get(&r.id).copied().unwrap_or(0.0);
        if coef < 0.0 {
            r.lb = -limit;
        } else {
            r.ub = limit;
        }
    }
    Ok(())
}

/// A minimal medium.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MinimalMedium {
    /// Exchange reaction and required uptake.
    pub uptake: BTreeMap<String, f64>,
    /// Total uptake, the objective of the linear form.
    pub total: f64,
}

/// Minimizes total uptake, or with `minimize_components` the number of imported exchanges first,
/// subject to the objective staying at `min_growth`. The integer form needs a backend with MILP.
pub fn minimal_medium(model: &Model, min_growth: f64, minimize_components: bool) -> Result<MinimalMedium> {
    let exchanges = model.exchanges();
    let n = model.reactions.len();
    let mut p: Problem = model.to_lp();
    p.sense = Some(Sense::Minimize);
    p.obj.iter_mut().for_each(|c| *c = 0.0);
    let coefs: Vec<(usize, f64)> = model.objective.clone();
    p.add_row(min_growth, f64::INFINITY, coefs);
    // Uptake variable u >= -(sign) v for every exchange.
    let mut uptake_cols = Vec::with_capacity(exchanges.len());
    for &j in &exchanges {
        let sign = model.exchange_coefficient(j).unwrap_or(-1.0).signum();
        let limit = if sign < 0.0 { -model.reactions[j].lb } else { model.reactions[j].ub };
        let u = p.add_col(1.0, 0.0, limit.max(0.0), false);
        // Uptake is flux in the direction that adds metabolite: u >= sign * v.
        p.add_row(0.0, f64::INFINITY, vec![(u, 1.0), (j, -sign)]);
        uptake_cols.push(u);
    }
    let backend = lp::default_backend();
    if minimize_components {
        let mut c = p.clone();
        c.obj.iter_mut().for_each(|o| *o = 0.0);
        let mut binaries = Vec::with_capacity(uptake_cols.len());
        for &u in &uptake_cols {
            let y = c.add_col(1.0, 0.0, 1.0, true);
            let limit = c.col_ub[u];
            let m = if limit.is_finite() { limit } else { 1e4 };
            c.add_row(f64::NEG_INFINITY, 0.0, vec![(u, 1.0), (y, -m)]);
            binaries.push(y);
        }
        let mut session = backend.open(&c)?;
        let s = session.solve()?;
        if !s.is_optimal() {
            return Err(Error::Solver(format!("minimal medium MILP ended with {:?}", s.status)));
        }
        // Fix the selection, then minimize total uptake over it.
        for (k, &y) in binaries.iter().enumerate() {
            if s.x[y] < 0.5 {
                p.col_ub[uptake_cols[k]] = 0.0;
            }
        }
    }
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    if !s.is_optimal() {
        return Err(Error::Solver(format!("minimal medium LP ended with {:?}", s.status)));
    }
    let mut uptake = BTreeMap::new();
    let mut total = 0.0;
    for (k, &j) in exchanges.iter().enumerate() {
        let value = s.x[uptake_cols[k]];
        total += value;
        if value > 1e-9 {
            uptake.insert(model.reactions[j].id.clone(), value);
        }
    }
    let _ = n;
    Ok(MinimalMedium { uptake, total })
}
