//! Dynamic FBA by static optimization (Mahadevan, Edwards and Doyle 2002; Varma and Palsson 1994).
//!
//! At each step the uptake limits follow Michaelis-Menten kinetics of the extracellular
//! concentrations, FBA gives the growth rate and exchange fluxes, and biomass and concentrations
//! advance over `dt` with those fluxes held constant:
//! `X(t+dt) = X exp(mu dt)`, `S_i(t+dt) = S_i + v_i X (exp(mu dt) - 1) / mu`. Uptake is capped so a
//! concentration cannot go below zero within the step.

use crate::error::{Error, Result};
use crate::lp::{self, Status};
use crate::model::Model;
use serde::{Deserialize, Serialize};

/// A substrate with Michaelis-Menten uptake.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Substrate {
    /// Exchange reaction identifier.
    pub exchange: String,
    /// Initial extracellular concentration (mmol/L).
    pub initial: f64,
    /// Maximum uptake (mmol/gDW/h).
    pub vmax: f64,
    /// Half-saturation constant (mmol/L).
    pub km: f64,
}

/// Settings of [`dfba`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DfbaOptions {
    /// Initial biomass (gDW/L).
    pub biomass: f64,
    /// End time (h).
    pub t_end: f64,
    /// Step (h).
    pub dt: f64,
    /// Substrates with kinetic uptake.
    pub substrates: Vec<Substrate>,
    /// Further exchange reactions whose concentrations are tracked, with initial values.
    pub tracked: Vec<(String, f64)>,
}

/// Time course of a dynamic FBA run.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DfbaResult {
    /// Sample times.
    pub time: Vec<f64>,
    /// Biomass at each time.
    pub biomass: Vec<f64>,
    /// Growth rate over the step starting at each time.
    pub growth_rate: Vec<f64>,
    /// Tracked exchange identifiers.
    pub species: Vec<String>,
    /// Concentrations: `concentrations[k][i]` is species `i` at time `k`.
    pub concentrations: Vec<Vec<f64>>,
}

/// Runs dynamic FBA.
pub fn dfba(model: &Model, options: &DfbaOptions) -> Result<DfbaResult> {
    if options.dt <= 0.0 || options.t_end <= 0.0 {
        return Err(Error::Invalid("dt and t_end must be positive".into()));
    }
    let backend = lp::default_backend();
    let mut session = backend.open(&model.to_lp())?;
    let mut ids = Vec::new();
    let mut cols = Vec::new();
    let mut conc = Vec::new();
    let mut kinetics = Vec::new();
    for s in &options.substrates {
        let j = model.reaction_index(&s.exchange)?;
        kinetics.push((ids.len(), j, s.vmax, s.km));
        ids.push(s.exchange.clone());
        cols.push(j);
        conc.push(s.initial);
    }
    for (id, initial) in &options.tracked {
        cols.push(model.reaction_index(id)?);
        ids.push(id.clone());
        conc.push(*initial);
    }
    let steps = (options.t_end / options.dt).round() as usize;
    let mut result = DfbaResult {
        time: Vec::with_capacity(steps + 1),
        biomass: Vec::with_capacity(steps + 1),
        growth_rate: Vec::with_capacity(steps + 1),
        species: ids,
        concentrations: Vec::with_capacity(steps + 1),
    };
    let mut x = options.biomass;
    for k in 0..=steps {
        let t = k as f64 * options.dt;
        for &(i, j, vmax, km) in &kinetics {
            let s = conc[i].max(0.0);
            let mut limit = vmax * s / (km + s);
            // Uptake that would exhaust the substrate within the step is reduced to what remains.
            limit = limit.min(s / (x.max(1e-300) * options.dt));
            session.set_col_bounds(j, -limit, model.reactions[j].ub)?;
        }
        let s = session.solve()?;
        let (mu, fluxes) = match s.status {
            Status::Optimal => (s.objective.max(0.0), Some(s.x)),
            _ => (0.0, None),
        };
        result.time.push(t);
        result.biomass.push(x);
        result.growth_rate.push(mu);
        result.concentrations.push(conc.clone());
        if k == steps {
            break;
        }
        if let Some(v) = fluxes {
            let integral = if mu > 1e-12 { x * ((mu * options.dt).exp() - 1.0) / mu } else { x * options.dt };
            for (i, &j) in cols.iter().enumerate() {
                conc[i] = (conc[i] + v[j] * integral).max(0.0);
            }
            x *= (mu * options.dt).exp();
        }
    }
    Ok(result)
}
