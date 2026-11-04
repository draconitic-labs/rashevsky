//! Dynamic FBA by static optimization (Mahadevan, Edwards and Doyle 2002; Varma and Palsson 1994).
//!
//! At each step the uptake limits follow Michaelis-Menten kinetics of the extracellular
//! concentrations, FBA gives the growth rate and exchange fluxes, and biomass and concentrations
//! advance over `dt` with those fluxes held constant:
//! `X(t+dt) = X exp(mu dt)`, `S_i(t+dt) = S_i - c_i v_i X (exp(mu dt) - 1) / mu`, with `c_i` the
//! exchange's coefficient on its metabolite, so that either way of writing an exchange works. The kinetic
//! limit never widens the model's own bounds, and what is there to take up overrides a bound that forces
//! uptake. The growth rate is the flux of the one reaction in the objective. Uptake is capped so that a
//! concentration stays above zero over the step: since the biomass integral grows with the growth rate
//! and the growth rate with the cap, the cap that uses up exactly what remains is found by bisection.

use crate::error::{Error, Result};
use crate::lp::{self, Sense, Session, Status};
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
            let lb = r.lb.max(-limit / per_flux);
            (lb, r.ub.max(lb))
        } else {
            let ub = r.ub.min(limit / per_flux);
            (r.lb.min(ub), ub)
        };
        session.set_col_bounds(self.column, lb, ub)
    }
}

/// The time integral of biomass over a step at growth rate `mu`.
fn biomass_integral(x: f64, mu: f64, dt: f64) -> f64 {
    if mu > 1e-12 { x * ((mu * dt).exp() - 1.0) / mu } else { x * dt }
}

/// What one solve of a step gives.
struct StepSolve {
    status: Status,
    mu: f64,
    fluxes: Option<Vec<f64>>,
}

/// Runs dynamic FBA.
pub fn dfba(model: &Model, options: &DfbaOptions) -> Result<DfbaResult> {
    if options.dt <= 0.0 || options.t_end <= 0.0 {
        return Err(Error::Invalid("dt and t_end must be positive".into()));
    }
    if !(options.biomass >= 0.0 && options.biomass.is_finite()) {
        return Err(Error::Invalid("the initial biomass must not be negative".into()));
    }
    let growth = match model.objective.as_slice() {
        [(j, coefficient)] if *coefficient > 0.0 && model.sense.unwrap_or(Sense::Maximize) == Sense::Maximize => *j,
        _ => {
            return Err(Error::Invalid(
                "dynamic FBA reads the growth rate from one maximized reaction, so the objective must be that".into(),
            ));
        }
    };
    let backend = lp::default_backend();
    let mut session = backend.open(&model.to_lp())?;
    let mut ids = Vec::new();
    let mut followed: Vec<Followed> = Vec::new();
    let mut conc = Vec::new();
    let mut follow = |id: &str, kinetics: Option<(f64, f64)>, initial: f64| -> Result<()> {
        let j = model.reaction_index(id)?;
        let coefficient =
            model.exchange_coefficient(j).ok_or_else(|| Error::Invalid(format!("{id} is not an exchange reaction")))?;
        if followed.iter().any(|f| f.column == j) {
            return Err(Error::Invalid(format!("{id} is followed twice")));
        }
        if !(initial >= 0.0 && initial.is_finite()) {
            return Err(Error::Invalid(format!("the concentration of {id} must not be negative")));
        }
        if kinetics.is_some_and(|(vmax, km)| vmax < 0.0 || km < 0.0) {
            return Err(Error::Invalid(format!("the kinetics of {id} must not be negative")));
        }
        followed.push(Followed { column: j, coefficient, kinetics });
        ids.push(id.to_string());
        conc.push(initial);
        Ok(())
    };
    for s in &options.substrates {
        follow(&s.exchange, Some((s.vmax, s.km)), s.initial)?;
    }
    for (id, initial) in &options.tracked {
        follow(id, None, *initial)?;
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
        let dt = if k == steps { options.dt } else { (options.t_end - t).min(options.dt) };
        let mut solve = |limits: &[f64]| -> Result<StepSolve> {
            for (f, &limit) in followed.iter().zip(limits) {
                f.apply(model, session.as_mut(), limit)?;
            }
            let s = session.solve()?;
            if s.status == Status::Unbounded {
                return Err(Error::Solver(format!(
                    "the flux problem is unbounded at t = {t}; give the exchange a bound"
                )));
            }
            Ok(if s.status == Status::Optimal {
                StepSolve { status: s.status, mu: s.x[growth].max(0.0), fluxes: Some(s.x) }
            } else {
                StepSolve { status: s.status, mu: 0.0, fluxes: None }
            })
        };
        // At zero growth the biomass integral is smallest, so this cap is the largest that can apply.
        let mut limits: Vec<f64> = followed
            .iter()
            .zip(&conc)
            .map(|(f, &s)| {
                let s = s.max(0.0);
                let kinetic = f.kinetics.map_or(f64::INFINITY, |(vmax, km)| vmax * s / (km + s).max(1e-300));
                kinetic.min(s / biomass_integral(x.max(1e-300), 0.0, dt))
            })
            .collect();
        let overdrawn = |i: usize, solved: &StepSolve, conc: &[f64]| -> bool {
            solved.fluxes.as_ref().is_some_and(|v| {
                let uptake = followed[i].coefficient * v[followed[i].column];
                uptake * biomass_integral(x, solved.mu, dt) > conc[i].max(0.0) * (1.0 + 1e-11) + 1e-13
            })
        };
        let mut solved = solve(&limits)?;
        for _ in 0..6 {
            let over: Vec<usize> = (0..followed.len()).filter(|&i| overdrawn(i, &solved, &conc)).collect();
            if over.is_empty() {
                break;
            }
            for i in over {
                let (mut low, mut high) = (0.0, limits[i]);
                for _ in 0..60 {
                    if high - low <= 1e-13 * high {
                        break;
                    }
                    limits[i] = 0.5 * (low + high);
                    let trial = solve(&limits)?;
                    if trial.fluxes.is_some() && !overdrawn(i, &trial, &conc) {
                        low = limits[i];
                    } else if trial.fluxes.is_none() {
                        // Too little to meet what the model requires, such as maintenance.
                        low = limits[i];
                    } else {
                        high = limits[i];
                    }
                }
                limits[i] = low;
            }
            solved = solve(&limits)?;
        }
        result.time.push(t);
        result.biomass.push(x);
        result.growth_rate.push(solved.mu);
        result.concentrations.push(conc.clone());
        result.status.push(solved.status);
        if k == steps {
            break;
        }
        if let Some(v) = &solved.fluxes {
            let integral = biomass_integral(x, solved.mu, dt);
            for (i, f) in followed.iter().enumerate() {
                conc[i] = (conc[i] - f.coefficient * v[f.column] * integral).max(0.0);
            }
            x *= (solved.mu * dt).exp();
        }
    }
    Ok(result)
}
