//! Dynamic FBA by static optimization (Mahadevan, Edwards and Doyle 2002; Varma and Palsson 1994).
//!
//! At each step the uptake limits follow Michaelis-Menten kinetics of the extracellular
//! concentrations, FBA gives the growth rate and exchange fluxes, and biomass and concentrations
//! advance over `dt` with those fluxes held constant:
//! `X(t+dt) = X exp(mu dt)`, `S_i(t+dt) = S_i - c_i v_i X (exp(mu dt) - 1) / mu`, with `c_i` the
//! exchange's coefficient on its metabolite, so that either way of writing an exchange works. Uptake is
//! capped by what the step's biomass integral can take without a concentration going below zero, and
//! the kinetic limit never widens the model's own bounds.

use crate::error::{Error, Result};
use crate::lp::{self, Session, Status};
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
    /// Solver status of the step starting at each time. Growth is zero when it is not optimal.
    pub status: Vec<Status>,
}

/// One exchange whose concentration is followed.
struct Followed {
    column: usize,
    /// Metabolite the exchange adds per unit of positive flux.
    coefficient: f64,
    /// `(vmax, km)` when uptake is kinetic.
    kinetics: Option<(f64, f64)>,
}

impl Followed {
    /// Sets flux bounds for an uptake of at most `limit`, within the model's own bounds.
    fn apply(&self, model: &Model, session: &mut dyn Session, limit: f64) -> Result<()> {
        let r = &model.reactions[self.column];
        let per_flux = self.coefficient.abs().max(1e-300);
        let (lb, ub) = if self.coefficient < 0.0 {
            (r.lb.max(-limit / per_flux), r.ub)
        } else {
            (r.lb, r.ub.min(limit / per_flux))
        };
        session.set_col_bounds(self.column, lb.min(ub), ub)
    }
}

/// The time integral of biomass over a step at growth rate `mu`.
fn biomass_integral(x: f64, mu: f64, dt: f64) -> f64 {
    if mu > 1e-12 { x * ((mu * dt).exp() - 1.0) / mu } else { x * dt }
}

/// Runs dynamic FBA.
pub fn dfba(model: &Model, options: &DfbaOptions) -> Result<DfbaResult> {
    if options.dt <= 0.0 || options.t_end <= 0.0 {
        return Err(Error::Invalid("dt and t_end must be positive".into()));
    }
    let backend = lp::default_backend();
    let mut session = backend.open(&model.to_lp())?;
    let mut ids = Vec::new();
    let mut followed = Vec::new();
    let mut conc = Vec::new();
    let coefficient = |j: usize, id: &str| {
        model.exchange_coefficient(j).ok_or_else(|| Error::Invalid(format!("{id} is not an exchange reaction")))
    };
    for s in &options.substrates {
        let j = model.reaction_index(&s.exchange)?;
        followed.push(Followed {
            column: j,
            coefficient: coefficient(j, &s.exchange)?,
            kinetics: Some((s.vmax, s.km)),
        });
        ids.push(s.exchange.clone());
        conc.push(s.initial);
    }
    for (id, initial) in &options.tracked {
        let j = model.reaction_index(id)?;
        followed.push(Followed { column: j, coefficient: coefficient(j, id)?, kinetics: None });
        ids.push(id.clone());
        conc.push(*initial);
    }
    let steps = ((options.t_end / options.dt) - 1e-9).ceil().max(1.0) as usize;
    let mut result = DfbaResult {
        time: Vec::with_capacity(steps + 1),
        biomass: Vec::with_capacity(steps + 1),
        growth_rate: Vec::with_capacity(steps + 1),
        species: ids,
        concentrations: Vec::with_capacity(steps + 1),
        status: Vec::with_capacity(steps + 1),
    };
    let mut x = options.biomass;
    for k in 0..=steps {
        let t = (k as f64 * options.dt).min(options.t_end);
        let dt = (options.t_end - t).min(options.dt).max(if k == steps { options.dt } else { 0.0 });
        // The cap depends on the growth rate through the biomass integral, and the growth rate on the
        // cap, so the two are iterated until no concentration would go below zero.
        let mut limits: Vec<f64> = followed
            .iter()
            .zip(&conc)
            .map(|(f, &s)| {
                let s = s.max(0.0);
                let kinetic = f.kinetics.map_or(f64::INFINITY, |(vmax, km)| vmax * s / (km + s));
                kinetic.min(s / biomass_integral(x.max(1e-300), 0.0, dt))
            })
            .collect();
        let mut status = Status::Unknown;
        let mut mu = 0.0;
        let mut fluxes = None;
        for _ in 0..40 {
            for (f, &limit) in followed.iter().zip(&limits) {
                f.apply(model, session.as_mut(), limit)?;
            }
            let s = session.solve()?;
            status = s.status;
            if status == Status::Unbounded {
                return Err(Error::Solver(format!(
                    "the flux problem is unbounded at t = {t}; give the exchange a bound"
                )));
            }
            if status != Status::Optimal {
                fluxes = None;
                mu = 0.0;
                break;
            }
            mu = s.objective.max(0.0);
            let integral = biomass_integral(x, mu, dt);
            let mut over = false;
            for (i, f) in followed.iter().enumerate() {
                let uptake = f.coefficient * s.x[f.column];
                let available = conc[i].max(0.0);
                if uptake * integral > available * (1.0 + 1e-12) + 1e-12 {
                    limits[i] = available / integral;
                    over = true;
                }
            }
            fluxes = Some(s.x);
            if !over {
                break;
            }
        }
        result.time.push(t);
        result.biomass.push(x);
        result.growth_rate.push(mu);
        result.concentrations.push(conc.clone());
        result.status.push(status);
        if k == steps {
            break;
        }
        if let Some(v) = fluxes {
            let integral = biomass_integral(x, mu, dt);
            for (i, f) in followed.iter().enumerate() {
                conc[i] = (conc[i] - f.coefficient * v[f.column] * integral).max(0.0);
            }
            x *= (mu * dt).exp();
        }
    }
    Ok(result)
}
