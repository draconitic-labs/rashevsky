//! Steady states of kinetic models with moiety conservation.
//!
//! The stoichiometric matrix `N` factors as `N = L Nr`, with `Nr` the independent rows and `L` the
//! link matrix (Reder 1988). Dependent species follow from the independent ones through the
//! conserved totals, and damped Newton solves `Nr v(y(s)) = 0` over the independent amounts `s`.
//!
//! Newton can reach any root from where integration stops, so a root that is unstable or has negative
//! amounts sends the integration further before it is accepted. The result reports the largest real
//! part of the Jacobian's eigenvalues: a positive one means the state is not an attractor, as on a
//! limit cycle. `delay` reads the present at a steady state, where eigenvalues of the delay-free system
//! do not decide stability, so a delayed model counts as stable only if the integration reached the root.
//! The integration extends while events still fire or are pending and while rates still change with the
//! time. If that does not end, there is no steady state, and that is an error, as are an event that would
//! fire on the way to the root and a Newton iteration that runs away.

use super::compile::{Coef, KineticModel, NoDelay, TimeSteps, eval};
use super::ode::OdeSystem;
use super::sim::{KineticSystem, Method, SimOptions, simulate};
use crate::error::{Error, Result};
use crate::linalg::rref;
use nalgebra::{DMatrix, DVector};

/// The factorization `N = L Nr` of the stoichiometric matrix over the state variables.
#[derive(Clone, Debug)]
pub struct Moieties {
    /// `states x reactions`.
    pub n: DMatrix<f64>,
    /// State indices of the independent rows.
    pub independent: Vec<usize>,
    /// `independent x reactions`.
    pub nr: DMatrix<f64>,
    /// `states x independent`.
    pub link: DMatrix<f64>,
}

/// Builds the stoichiometric matrix over the state variables with coefficients at their initial values.
pub fn stoichiometric_matrix(model: &KineticModel) -> DMatrix<f64> {
    let mut n = DMatrix::zeros(model.n_states(), model.reactions.len());
    for e in &model.stoich {
        let coef = match e.coef {
            Coef::Const(v) => v,
            Coef::Slot(s) => model.initial[s],
        };
        let conv = e.conversion.map_or(1.0, |s| model.initial[s]);
        n[(e.state, e.reaction)] += e.sign * coef * conv;
    }
    n
}

/// Finds independent rows and the link matrix.
pub fn moieties(model: &KineticModel) -> Result<Moieties> {
    if model.has_varying_stoichiometry() {
        return Err(Error::Kinetic(
            "moiety analysis needs fixed stoichiometry, and a coefficient or conversion factor here varies".into(),
        ));
    }
    if !model.rate_rules.is_empty() || !model.algebraic.is_empty() {
        return Err(Error::Kinetic(
            "moiety analysis applies to reaction networks without rate or algebraic rules".into(),
        ));
    }
    let n = stoichiometric_matrix(model);
    let e = rref(&n.transpose(), 1e-10);
    let mut independent = e.pivot_cols.clone();
    independent.sort_unstable();
    let r = independent.len();
    let mut nr = DMatrix::zeros(r, n.ncols());
    for (i, &row) in independent.iter().enumerate() {
        nr.row_mut(i).copy_from(&n.row(row));
    }
    let gram = &nr * nr.transpose();
    let inv = gram.lu().try_inverse().ok_or_else(|| Error::Kinetic("singular independent rows".into()))?;
    let link = &n * nr.transpose() * inv;
    Ok(Moieties { n, independent, nr, link })
}

/// Settings of [`steady_state`].
#[derive(Clone, Copy, Debug)]
pub struct SteadyOptions {
    /// Largest entry of the rate of change accepted at the solution.
    pub tolerance: f64,
    /// Integrate over this time first, which selects the attractor reached from the initial state. It
    /// is extended past the last time a rate can jump.
    pub integrate_first: f64,
    /// Newton iterations.
    pub max_iterations: usize,
}

impl Default for SteadyOptions {
    fn default() -> Self {
        Self { tolerance: 1e-10, integrate_first: 100.0, max_iterations: 200 }
    }
}

/// A steady state.
#[derive(Clone, Debug)]
pub struct SteadyState {
    /// State at the steady state.
    pub y: Vec<f64>,
    /// Slot values at the steady state.
    pub values: Vec<f64>,
    /// Flux per reaction.
    pub fluxes: Vec<f64>,
    /// Largest absolute rate of change.
    pub residual: f64,
    /// Newton iterations used.
    pub iterations: usize,
    /// Time at which the rates were evaluated.
    pub time: f64,
    /// Largest real part of the eigenvalues of the Jacobian over the independent amounts. NaN for a model
    /// with `delay`, whose spectrum this does not describe.
    pub largest_eigenvalue: f64,
    /// Whether every perturbation decays. For a model with `delay`, whether the integration reached this
    /// state.
    pub stable: bool,
}

/// How many times the integration is extended, tenfold each, in search of an attractor. A limit cycle uses all of them.
const EXTENSIONS: usize = 2;

/// Damped Newton from `start`, evaluating at time `time`.
fn newton(
    model: &KineticModel,
    m: &Moieties,
    start: &[f64],
    time: f64,
    options: &SteadyOptions,
) -> Result<SteadyState> {
    let sys = KineticSystem::new(model, false);
    let n = start.len();
    let y0 = start.to_vec();
    let ind = &m.independent;
    let lift = |s: &[f64]| -> Vec<f64> {
        // y = y0 + L (s - s0)
        let ds = DVector::from_iterator(s.len(), s.iter().zip(ind).map(|(a, &i)| a - y0[i]));
        let dy = &m.link * ds;
        (0..n).map(|k| y0[k] + dy[k]).collect()
    };
    let residual_of = |s: &[f64]| -> (Vec<f64>, f64) {
        let full = lift(s);
        let mut f = vec![0.0; n];
        sys.rhs(time, &full, &mut f);
        let norm = f.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
        (f, norm)
    };
    let reduced_jacobian = |s: &[f64]| -> DMatrix<f64> {
        let full = lift(s);
        let mut jac = DMatrix::zeros(n, n);
        sys.jacobian(time, &full, &mut jac);
        let mut reduced = DMatrix::zeros(ind.len(), ind.len());
        let jl = &jac * &m.link;
        for (row, &i) in ind.iter().enumerate() {
            reduced.row_mut(row).copy_from(&jl.row(i));
        }
        reduced
    };
    let mut s: Vec<f64> = ind.iter().map(|&i| y0[i]).collect();
    let (mut f, mut norm) = residual_of(&s);
    let mut iterations = 0;
    // Newton runs until it stops improving, not merely to the tolerance: control coefficients divide
    // by small fluxes, so they need the steady state as exact as the arithmetic allows.
    while norm > 0.0 && iterations < options.max_iterations {
        iterations += 1;
        let reduced = reduced_jacobian(&s);
        let rhs = DVector::from_iterator(ind.len(), ind.iter().map(|&i| -f[i]));
        let Some(step) = reduced.lu().solve(&rhs) else {
            return Err(Error::Kinetic("singular Jacobian in the steady-state iteration".into()));
        };
        let mut scale = 1.0;
        let mut accepted = false;
        for _ in 0..30 {
            let trial: Vec<f64> = s.iter().zip(step.iter()).map(|(a, d)| a + scale * d).collect();
            let (ft, nt) = residual_of(&trial);
            if nt.is_finite() && nt < norm {
                s = trial;
                f = ft;
                norm = nt;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    if norm > options.tolerance {
        return Err(Error::Kinetic(format!("no steady state found; residual {norm:e}")));
    }
    let full = lift(&s);
    let values = sys.values(time, &full);
    let fluxes = model.reactions.iter().map(|r| values[r.rate]).collect();
    // A root has a last Newton step that is small. Where the rate only tends to zero as an amount grows
    // without limit, as with `1 / (1 + S)`, the residual falls but the step does not.
    let (largest_eigenvalue, size) = if ind.is_empty() {
        (0.0, 0.0)
    } else {
        let reduced = reduced_jacobian(&s);
        if norm > 0.0 {
            let rhs = DVector::from_iterator(ind.len(), ind.iter().map(|&i| -f[i]));
            if let Some(step) = reduced.clone().lu().solve(&rhs) {
                let largest_step = step.iter().fold(0.0f64, |a, v| a.max(v.abs()));
                let largest_amount = s.iter().fold(0.0f64, |a, v| a.max(v.abs()));
                if largest_step > 1e-6 * (1.0 + largest_amount) {
                    return Err(Error::Kinetic(format!(
                        "no steady state found: the iteration runs away, with amounts near {largest_amount:e} and steps of {largest_step:e}"
                    )));
                }
            }
        }
        let eigenvalues = reduced.complex_eigenvalues();
        let largest = eigenvalues.iter().map(|z| z.re).fold(f64::NEG_INFINITY, f64::max);
        (largest, eigenvalues.iter().map(|z| z.norm()).fold(0.0f64, f64::max))
    };
    Ok(SteadyState {
        y: full,
        values,
        fluxes,
        residual: norm,
        iterations,
        time,
        largest_eigenvalue,
        stable: largest_eigenvalue <= 1e-9 * size,
    })
}

/// True when a species at this state has a clearly negative amount.
fn has_negative_amount(model: &KineticModel, state: &SteadyState) -> bool {
    let scale = model.species.iter().map(|s| state.values[s.amount].abs()).fold(0.0f64, f64::max).max(1.0);
    model.species.iter().any(|s| state.values[s.amount] < -1e-9 * scale)
}

/// Rates of change and reaction rates at a time, for a frozen state.
fn rates_at(model: &KineticModel, sys: &KineticSystem<'_>, y: &[f64], t: f64) -> Vec<f64> {
    let mut dy = vec![0.0; y.len()];
    sys.rhs(t, y, &mut dy);
    let vals = sys.values(t, y);
    dy.extend(model.reactions.iter().map(|r| vals[r.rate]));
    dy
}

/// Whether the rates at `y` no longer change with the time after `t`: no jump lies ahead on a long grid, and
/// they read the same at times spread over decades, compared relative to their own size.
fn time_has_settled(model: &KineticModel, steps: &TimeSteps, y: &[f64], t: f64) -> Result<bool> {
    if !model.reads_time() {
        return Ok(true);
    }
    let sys = KineticSystem::new(model, false);
    let vals = sys.values(t, y);
    let far = 100.0 * t + 1000.0;
    if !steps.rate_times(&vals, t, far)?.is_empty() || !steps.trigger_times(&vals, t, far)?.is_empty() {
        return Ok(false);
    }
    let now = rates_at(model, &sys, y, t);
    let mut offset = 0.37;
    while offset < 1e7 + t * 1e4 {
        let later = rates_at(model, &sys, y, t + offset);
        let size = now.iter().chain(&later).fold(0.0f64, |acc, v| acc.max(v.abs()));
        if now.iter().zip(&later).any(|(a, b)| (a - b).abs() > 1e-8 * size) {
            return Ok(false);
        }
        offset *= 1.73;
    }
    Ok(true)
}

/// Rejects a root that an event would stop the system from reaching.
fn check_events_quiet(analysed: &KineticModel, integrated: &[f64], root: &SteadyState) -> Result<()> {
    let model = analysed;
    if model.events.is_empty() {
        return Ok(());
    }
    let mut vals = model.initial.clone();
    model.eval_state(root.time, integrated, &mut vals, &NoDelay);
    for e in &model.events {
        let holds = |v: &[f64]| {
            let x = eval(&e.trigger, v, root.time, &NoDelay);
            x != 0.0 && !x.is_nan()
        };
        if holds(&vals) != holds(&root.values) {
            return Err(Error::Kinetic(format!(
                "event {} fires between the integrated state and the root, so that root is not reached",
                e.id
            )));
        }
    }
    Ok(())
}

/// Solves for a steady state.
pub fn steady_state(model: &KineticModel, options: &SteadyOptions) -> Result<SteadyState> {
    let m = moieties(model)?;
    let frozen = model.has_delay.then(|| model.without_delays());
    let analysed = frozen.as_ref().unwrap_or(model);
    let steps = analysed.time_steps();
    let integrate = |until: f64| -> Result<(Vec<f64>, Vec<f64>, Option<f64>)> {
        if until <= 0.0 {
            return Ok((model.initial_state(), Vec::new(), None));
        }
        let run = simulate(
            model,
            &[0.0, until],
            &SimOptions { method: Method::Rodas4, rtol: 1e-8, atol: 1e-10, ..SimOptions::default() },
        )?;
        let vals = &run.values[1];
        Ok((model.state_slots.iter().map(|&s| vals[s]).collect(), run.event_times, run.pending_until))
    };
    // Everything that happens at fixed times must have happened by the end of the integration.
    let reach = options.integrate_first.max(1.0);
    let last_jump = steps
        .rate_times(&model.initial, 0.0, reach)?
        .into_iter()
        .chain(steps.trigger_times(&model.initial, 0.0, reach)?)
        .fold(0.0f64, f64::max);
    let mut horizon = if options.integrate_first > 0.0 { options.integrate_first.max(2.0 * last_jump) } else { 0.0 };
    let mut fallback = None;
    let mut last_error = None;
    let mut extensions = 0;
    let mut extend = |horizon: &mut f64, reason: String| -> Result<()> {
        if *horizon <= 0.0 || extensions == EXTENSIONS {
            return Err(Error::Kinetic(reason));
        }
        extensions += 1;
        *horizon *= 10.0;
        Ok(())
    };
    for _ in 0..12 {
        // A run that diverges during an extension has told us all it can: the root found before is unstable.
        let (y, event_times, pending) = match integrate(horizon) {
            Ok(run) => run,
            Err(e) => return fallback.ok_or(e),
        };
        if let Some(until) = pending {
            horizon = horizon.max(1.5 * until + 1.0);
            continue;
        }
        if event_times.iter().any(|&t| t > 0.5 * horizon) {
            let reason = format!("events are still firing at t = {horizon}, so there is no steady state");
            extend(&mut horizon, reason)?;
            continue;
        }
        if !time_has_settled(analysed, &steps, &y, horizon)? {
            let reason = format!("the rates still depend on the time at t = {horizon}, so there is no steady state");
            extend(&mut horizon, reason)?;
            continue;
        }
        match newton(analysed, &m, &y, horizon, options) {
            Ok(mut state) => {
                check_events_quiet(analysed, &y, &state)?;
                if model.has_delay {
                    let distance = y.iter().zip(&state.y).fold(0.0f64, |a, (p, q)| a.max((p - q).abs()));
                    let size = state.y.iter().fold(0.0f64, |a, v| a.max(v.abs()));
                    state.stable = distance <= 1e-4 * size + 1e-10;
                    state.largest_eigenvalue = f64::NAN;
                }
                let negative = has_negative_amount(model, &state);
                if state.stable && !negative {
                    return Ok(state);
                }
                if !negative {
                    fallback = Some(state);
                }
            }
            Err(e) => {
                if horizon <= 0.0 {
                    return Err(e);
                }
                last_error = Some(e);
            }
        }
        if extend(&mut horizon, String::new()).is_err() {
            break;
        }
    }
    if let Some(e) = last_error
        && fallback.is_none()
    {
        return Err(e);
    }
    fallback.ok_or_else(|| Error::Kinetic("the only steady state found has negative amounts".into()))
}
