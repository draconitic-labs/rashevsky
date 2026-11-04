//! Time-course simulation with SBML event semantics.
//!
//! Between events the system integrates with a stiff method. Trigger transitions are bracketed on
//! each step, located by bisection on the step's interpolant, and the step is then repeated exactly
//! to the located time. Event execution follows SBML Level 3: triggers fire on a false-to-true
//! transition, delays and priorities are evaluated as specified, each executed event re-evaluates all
//! triggers, and integration restarts from the new state.

use super::compile::{C, DelayContext, KineticModel, PreStart, TimeSteps, eval};
use super::ode::{Bdf, Integrator, OdeSystem, Rodas4, Stats, Tolerances};
use crate::error::{Error, Result};
use nalgebra::DMatrix;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;

/// Integration method.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Method {
    /// Rosenbrock order 4(3).
    #[default]
    Rodas4,
    /// Variable-order backward differences.
    Bdf,
}

/// Settings of a simulation.
#[derive(Clone, Copy, Debug)]
pub struct SimOptions {
    /// Integration method.
    pub method: Method,
    /// Relative tolerance.
    pub rtol: f64,
    /// Absolute tolerance.
    pub atol: f64,
    /// Use the numeric Jacobian even when the analytic one exists.
    pub numeric_jacobian: bool,
    /// Seed of the generator that breaks ties between events of equal priority.
    pub seed: u64,
}

impl Default for SimOptions {
    fn default() -> Self {
        Self { method: Method::Rodas4, rtol: 1e-9, atol: 1e-12, numeric_jacobian: false, seed: 1 }
    }
}

/// Simulated values of every symbol at the requested times.
#[derive(Clone, Debug)]
pub struct TimeCourse {
    /// Output times.
    pub time: Vec<f64>,
    /// `values[k][slot]` is the value of a slot at `time[k]`.
    pub values: Vec<Vec<f64>>,
    /// Integrator work.
    pub stats: Stats,
    /// When each event execution happened.
    pub event_times: Vec<f64>,
    /// When the latest event still pending at the end of the run would execute.
    pub pending_until: Option<f64>,
}

impl TimeCourse {
    /// The series of a named symbol: a slot name, or `S#a` for the amount of species `S`.
    pub fn series(&self, model: &KineticModel, name: &str) -> Result<Vec<f64>> {
        let slot = model.slot(name)?;
        Ok(self.values.iter().map(|v| v[slot]).collect())
    }

    /// The concentration series of a species.
    pub fn concentration(&self, model: &KineticModel, species: &str) -> Result<Vec<f64>> {
        let info = model
            .species
            .iter()
            .find(|s| s.id == species)
            .ok_or_else(|| Error::Unknown { kind: "species", id: species.into() })?;
        Ok(self
            .values
            .iter()
            .map(|v| if info.has_only_substance_units { v[info.amount] / v[info.compartment] } else { v[info.view] })
            .collect())
    }
}

#[derive(Clone)]
struct HistoryPoint {
    t: f64,
    y: Vec<f64>,
    f: Vec<f64>,
    vals: Vec<f64>,
}

/// The kinetic model as an ODE system, with the solution history that `delay` reads.
pub struct KineticSystem<'a> {
    model: &'a KineticModel,
    vals: RefCell<Vec<f64>>,
    history: RefCell<Vec<HistoryPoint>>,
    t0: f64,
    y0: Vec<f64>,
    vals0: Vec<f64>,
    numeric_jacobian: bool,
}

impl<'a> KineticSystem<'a> {
    /// Wraps a model.
    pub fn new(model: &'a KineticModel, numeric_jacobian: bool) -> Self {
        Self {
            model,
            vals: RefCell::new(model.initial.clone()),
            history: RefCell::new(Vec::new()),
            t0: 0.0,
            y0: model.initial_state(),
            vals0: model.initial.clone(),
            // The symbolic derivatives leave `delay` out, which would cost the integrator its order.
            numeric_jacobian: numeric_jacobian || model.jacobian.is_none() || model.has_delay,
        }
    }

    fn interpolated_state(&self, t: f64) -> Vec<f64> {
        let h = self.history.borrow();
        if t <= self.t0 || h.is_empty() {
            return self.y0.clone();
        }
        // The first point with time >= t ends the bracketing interval.
        let idx = h.partition_point(|p| p.t < t);
        if idx == 0 {
            return h[0].y.clone();
        }
        if idx >= h.len() {
            return h[h.len() - 1].y.clone();
        }
        let (a, b) = (&h[idx - 1], &h[idx]);
        let step = b.t - a.t;
        if step <= 0.0 {
            return b.y.clone();
        }
        let s = (t - a.t) / step;
        let (h00, h10) = ((1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s), s * (1.0 - s) * (1.0 - s));
        let (h01, h11) = (s * s * (3.0 - 2.0 * s), s * s * (s - 1.0));
        (0..a.y.len()).map(|k| h00 * a.y[k] + h10 * step * a.f[k] + h01 * b.y[k] + h11 * step * b.f[k]).collect()
    }

    fn pre_start(&self) -> PreStart<'_> {
        PreStart {
            items: &self.model.initial_items,
            base: &self.model.initial,
            rate_of_slots: &self.model.rate_of_slots,
        }
    }

    /// All slots at a past time, by linear interpolation between stored steps.
    fn stored_values(&self, t: f64) -> Vec<f64> {
        let h = self.history.borrow();
        if t <= self.t0 || h.is_empty() {
            return self.vals0.clone();
        }
        let idx = h.partition_point(|p| p.t < t);
        if idx == 0 {
            return h[0].vals.clone();
        }
        if idx >= h.len() {
            return h[h.len() - 1].vals.clone();
        }
        let (a, b) = (&h[idx - 1], &h[idx]);
        let w = if b.t > a.t { (t - a.t) / (b.t - a.t) } else { 1.0 };
        a.vals.iter().zip(&b.vals).map(|(x, y)| x + w * (y - x)).collect()
    }

    fn fill_vals(&self, t: f64, y: &[f64]) -> Vec<f64> {
        let mut vals = std::mem::take(&mut *self.vals.borrow_mut());
        self.model.eval_state(t, y, &mut vals, self);
        vals
    }

    fn put_vals(&self, vals: Vec<f64>) {
        *self.vals.borrow_mut() = vals;
    }

    /// All slot values for a state.
    pub fn values(&self, t: f64, y: &[f64]) -> Vec<f64> {
        let mut vals = self.model.initial.clone();
        self.model.eval_state(t, y, &mut vals, self);
        vals
    }
}

/// Reads `delay` inside an already delayed expression from stored values, so that nested delays
/// cost one interpolation each rather than a recursive re-evaluation of the model.
struct StoredDelay<'s, 'a> {
    sys: &'s KineticSystem<'a>,
}

impl DelayContext for StoredDelay<'_, '_> {
    fn value_at(&self, expr: &C, t: f64) -> f64 {
        if t < self.sys.t0 {
            return self.sys.pre_start().value_at(expr, t);
        }
        let vals = self.sys.stored_values(t);
        eval(expr, &vals, t, self)
    }
}

impl DelayContext for KineticSystem<'_> {
    fn value_at(&self, expr: &C, t: f64) -> f64 {
        if t < self.t0 {
            return self.pre_start().value_at(expr, t);
        }
        let y = self.interpolated_state(t);
        let mut vals = self.model.initial.clone();
        let shallow = StoredDelay { sys: self };
        self.model.eval_state(t, &y, &mut vals, &shallow);
        eval(expr, &vals, t, &shallow)
    }
}

impl OdeSystem for KineticSystem<'_> {
    fn dim(&self) -> usize {
        self.model.n_states()
    }

    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
        let vals = self.fill_vals(t, y);
        self.model.rhs_from_vals(t, &vals, dy, self);
        self.put_vals(vals);
    }

    fn jacobian(&self, t: f64, y: &[f64], jac: &mut DMatrix<f64>) {
        jac.fill(0.0);
        match (&self.model.jacobian, self.numeric_jacobian) {
            (Some(entries), false) => {
                let vals = self.fill_vals(t, y);
                for (i, j, c) in entries {
                    jac[(*i, *j)] += eval(c, &vals, t, self);
                }
                self.put_vals(vals);
            }
            _ => {
                let n = y.len();
                let mut base = vec![0.0; n];
                self.rhs(t, y, &mut base);
                let mut plus = vec![0.0; n];
                let mut probe = y.to_vec();
                for j in 0..n {
                    let h = (f64::EPSILON * y[j].abs().max(1e-5)).sqrt();
                    probe[j] = y[j] + h;
                    self.rhs(t, &probe, &mut plus);
                    probe[j] = y[j];
                    for i in 0..n {
                        jac[(i, j)] = (plus[i] - base[i]) / h;
                    }
                }
            }
        }
    }

    fn time_derivative(&self, t: f64, y: &[f64], out: &mut [f64]) {
        out.fill(0.0);
        if self.numeric_jacobian {
            let n = y.len();
            let mut a = vec![0.0; n];
            let mut b = vec![0.0; n];
            let h = (f64::EPSILON * t.abs().max(1e-5)).sqrt();
            self.rhs(t, y, &mut a);
            self.rhs(t + h, y, &mut b);
            for i in 0..n {
                out[i] = (b[i] - a[i]) / h;
            }
            return;
        }
        let vals = self.fill_vals(t, y);
        for (i, c) in &self.model.time_derivative {
            out[*i] = eval(c, &vals, t, self);
        }
        self.put_vals(vals);
    }

    fn mass(&self) -> &[f64] {
        &self.model.mass
    }
}

/// Points per step at which triggers are examined for a transition.
const TRIGGER_SCAN: usize = 16;

#[derive(Clone, Debug)]
struct PendingEvent {
    event: usize,
    exec_time: f64,
    values: Option<Vec<f64>>,
}

struct Driver<'a> {
    model: &'a KineticModel,
    sys: KineticSystem<'a>,
    int: Box<dyn Integrator>,
    tol: Tolerances,
    trig: Vec<bool>,
    pending: Vec<PendingEvent>,
    rng: ChaCha8Rng,
    /// Times at which a rate can jump, and at which a trigger can change.
    rate_thresholds: Vec<f64>,
    trigger_thresholds: Vec<f64>,
    steps: TimeSteps,
    /// Last output time.
    end: f64,
    event_times: Vec<f64>,
}

/// How far either side of a rate threshold integration stops. The state moves by a negligible amount
/// over this gap, and no stage of a step sees the jump.
fn threshold_margin(t: f64) -> f64 {
    4e-14 * (1.0 + t.abs())
}

impl Driver<'_> {
    fn triggers_at(&self, t: f64, y: &[f64]) -> Vec<bool> {
        let vals = self.sys.values(t, y);
        self.model
            .events
            .iter()
            .map(|e| {
                let v = eval(&e.trigger, &vals, t, &self.sys);
                v != 0.0 && !v.is_nan()
            })
            .collect()
    }

    /// True when the trigger vector differs from the stored one in a way that needs action.
    fn transition(&self, current: &[bool]) -> bool {
        current.iter().zip(&self.trig).any(|(a, b)| a != b)
    }

    /// Fires transitions and executes everything due at time `t`. Returns the state after.
    fn handle_events(&mut self, t: f64, y: &mut [f64]) -> Result<bool> {
        let mut changed = false;
        let tiny = 1e-12 * (1.0 + t.abs());
        for _guard in 0..100_000 {
            let vals = self.sys.values(t, y);
            let current = self.triggers_at(t, y);
            for (k, e) in self.model.events.iter().enumerate() {
                if current[k] && !self.trig[k] {
                    let delay = e.delay.as_ref().map_or(0.0, |d| eval(d, &vals, t, &self.sys));
                    if delay < 0.0 || delay.is_nan() {
                        return Err(Error::Kinetic(format!("event {} has delay {delay} at t = {t}", e.id)));
                    }
                    let values = if e.use_values_from_trigger_time {
                        Some(e.assignments.iter().map(|(_, c)| eval(c, &vals, t, &self.sys)).collect())
                    } else {
                        None
                    };
                    self.pending.push(PendingEvent { event: k, exec_time: t + delay, values });
                } else if !current[k] && self.trig[k] && !e.persistent {
                    self.pending.retain(|p| p.event != k);
                }
            }
            self.trig = current;
            // Pick the highest-priority event due now.
            // Ties are broken uniformly at random, as SBML specifies, from a seeded generator.
            let mut best: Option<(usize, f64)> = None;
            let mut ties = 1u32;
            for (i, p) in self.pending.iter().enumerate() {
                if p.exec_time <= t + tiny {
                    let e = &self.model.events[p.event];
                    let priority = e.priority.as_ref().map_or(f64::NEG_INFINITY, |c| eval(c, &vals, t, &self.sys));
                    match best {
                        None => best = Some((i, priority)),
                        Some((_, bp)) if priority > bp => {
                            best = Some((i, priority));
                            ties = 1;
                        }
                        Some((_, bp)) if priority == bp => {
                            ties += 1;
                            if self.rng.random_range(0..ties) == 0 {
                                best = Some((i, priority));
                            }
                        }
                        _ => {}
                    }
                }
            }
            let Some((i, _)) = best else {
                if changed && self.steps.dynamic {
                    self.find_thresholds(t, y)?;
                }
                return Ok(changed);
            };
            let due = self.pending.remove(i);
            self.event_times.push(t);
            let e = &self.model.events[due.event];
            let values: Vec<f64> = match due.values {
                Some(v) => v,
                None => e.assignments.iter().map(|(_, c)| eval(c, &vals, t, &self.sys)).collect(),
            };
            self.model.apply_event(due.event, &values, t, y, &self.sys);
            self.project_algebraic(t, y)?;
            changed = true;
        }
        Err(Error::Kinetic(format!("event loop did not settle at t = {t}")))
    }

    /// Restores the algebraic rules after a jump, by Newton iteration on the variables they determine.
    fn project_algebraic(&self, t: f64, y: &mut [f64]) -> Result<()> {
        let m = self.model;
        if m.algebraic.is_empty() {
            return Ok(());
        }
        let n = m.algebraic.len();
        let residual = |y: &[f64]| -> Vec<f64> {
            let vals = self.sys.values(t, y);
            m.algebraic.iter().map(|(_, c)| eval(c, &vals, t, &self.sys)).collect()
        };
        for _ in 0..20 {
            let r = residual(y);
            let norm = r.iter().map(|v| v * v).sum::<f64>().sqrt();
            if norm < 1e-13 {
                return Ok(());
            }
            if !norm.is_finite() {
                break;
            }
            let mut jac = DMatrix::zeros(n, n);
            for (col, (state, _)) in m.algebraic.iter().enumerate() {
                let save = y[*state];
                let h = 1e-7 * save.abs().max(1e-3);
                y[*state] = save + h;
                let plus = residual(y);
                y[*state] = save;
                for row in 0..n {
                    jac[(row, col)] = (plus[row] - r[row]) / h;
                }
            }
            let rhs = nalgebra::DVector::from_iterator(n, r.iter().map(|v| -v));
            let Some(dx) = jac.lu().solve(&rhs) else { break };
            for (k, (state, _)) in m.algebraic.iter().enumerate() {
                y[*state] += dx[k];
            }
        }
        let remaining = residual(y).iter().map(|v| v * v).sum::<f64>().sqrt();
        if remaining < 1e-8 * (1.0 + y.iter().fold(0.0f64, |a, v| a.max(v.abs()))) {
            Ok(())
        } else {
            Err(Error::Kinetic(format!(
                "the algebraic rules cannot be satisfied after the event at t = {t}; residual {remaining:e}"
            )))
        }
    }

    /// Brings fast reactions to equilibrium with the time frozen.
    fn equilibrate(&mut self, t: f64, y: &mut [f64]) -> Result<()> {
        if !self.model.has_fast() {
            return Ok(());
        }
        let fast = FastOnly { sys: &self.sys, t, mass: vec![1.0; y.len()] };
        let tol = Tolerances { rtol: 1e-10, atol: 1e-13, ..Tolerances::default() };
        let mut int = Rodas4::new();
        int.restart(&fast, 0.0, y);
        let horizon = 1e-2;
        while int.time() < horizon {
            int.step(&fast, horizon, &tol)?;
        }
        y.copy_from_slice(int.state());
        Ok(())
    }

    fn next_exec(&self) -> f64 {
        self.pending.iter().map(|p| p.exec_time).fold(f64::INFINITY, f64::min)
    }

    /// Computes the times of jumps from the current values.
    fn find_thresholds(&mut self, t: f64, y: &[f64]) -> Result<()> {
        let vals = self.sys.values(t, y);
        self.rate_thresholds = self.steps.rate_times(&vals, t, self.end)?;
        self.trigger_thresholds = self.steps.trigger_times(&vals, t, self.end)?;
        Ok(())
    }

    /// Integrates to `limit`, keeping the history that `delay` reads.
    fn integrate_to(&mut self, limit: f64) -> Result<()> {
        while self.int.time() < limit {
            self.int.step(&self.sys, limit, &self.tol)?;
            if self.model.has_delay {
                let (t, y) = (self.int.time(), self.int.state().to_vec());
                self.push_history(t, &y);
            }
        }
        Ok(())
    }

    /// Moves to a time beside a rate threshold, over a gap small enough that the state does not change,
    /// and restarts integration there so that no stage of a step straddles the jump.
    fn move_to(&mut self, at: f64) -> Result<()> {
        let mut y = self.int.state().to_vec();
        self.push_history(at, &y);
        // An algebraic variable can jump with the rate.
        self.project_algebraic(at, &mut y)?;
        if self.handle_events(at, &mut y)? {
            self.equilibrate(at, &mut y)?;
            self.push_history(at, &y);
        }
        self.int.restart(&self.sys, at, &y);
        Ok(())
    }

    /// Integrates to `target`, handling every event on the way.
    fn advance(&mut self, target: f64) -> Result<()> {
        while self.int.time() < target {
            let t0 = self.int.time();
            let y0 = self.int.state().to_vec();
            let mut limit = target.min(self.next_exec()).max(t0);
            let next = self.rate_thresholds.partition_point(|&th| th + threshold_margin(th) <= t0);
            if let Some(&threshold) = self.rate_thresholds.get(next) {
                let below = threshold - threshold_margin(threshold);
                if t0 >= below {
                    if t0 < threshold {
                        // Reach the threshold itself without a step, whose last stage would see the jump.
                        self.move_to(threshold)?;
                    } else if target > threshold {
                        self.move_to(threshold + threshold_margin(threshold))?;
                    } else {
                        break;
                    }
                    continue;
                }
                limit = limit.min(below);
            }
            let next = self.trigger_thresholds.partition_point(|&th| th <= t0);
            if let Some(&threshold) = self.trigger_thresholds.get(next) {
                limit = limit.min(threshold);
            }
            if limit > t0 {
                self.int.step(&self.sys, limit, &self.tol)?;
            }
            let t1 = self.int.time();
            let y1 = self.int.state().to_vec();
            if self.model.has_delay {
                self.push_history(t1, &y1);
            }
            // Scan the step: a trigger window can open and close between its two ends.
            let mut bracket: Option<(f64, f64)> = None;
            let mut probe = vec![0.0; y1.len()];
            if !self.model.events.is_empty() {
                let mut previous = t0;
                for k in 1..=TRIGGER_SCAN {
                    let at = t0 + (t1 - t0) * k as f64 / TRIGGER_SCAN as f64;
                    let scan = if k == TRIGGER_SCAN {
                        self.triggers_at(t1, &y1)
                    } else {
                        self.int.interpolate(at, &mut probe);
                        self.triggers_at(at, &probe)
                    };
                    if self.transition(&scan) {
                        bracket = Some((previous, at));
                        break;
                    }
                    previous = at;
                }
            }
            if let Some((mut lo, mut hi)) = bracket {
                // Refine the first transition with the step's interpolant.
                while hi - lo > 1e-11 * (1.0 + hi.abs()) {
                    let mid = 0.5 * (lo + hi);
                    self.int.interpolate(mid, &mut probe);
                    let at_mid = self.triggers_at(mid, &probe);
                    if self.transition(&at_mid) {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                }
                // Repeat the step exactly up to the located time.
                self.int.restart(&self.sys, t0, &y0);
                if self.model.has_delay {
                    self.sys.history.borrow_mut().retain(|p| p.t <= t0);
                }
                self.integrate_to(hi)?;
                let mut y = self.int.state().to_vec();
                let located = self.triggers_at(hi, &y);
                if !self.transition(&located) {
                    // The exact state lags the interpolant; nudge across the transition.
                    let nudge = hi + 1e-10 * (1.0 + hi.abs());
                    if nudge < t1 {
                        self.integrate_to(nudge)?;
                        y = self.int.state().to_vec();
                    } else {
                        self.int.restart(&self.sys, t0, &y0);
                        if self.model.has_delay {
                            self.sys.history.borrow_mut().retain(|p| p.t <= t0);
                        }
                        self.integrate_to(t1)?;
                        y = y1.clone();
                    }
                }
                let at = self.int.time();
                self.handle_events(at, &mut y)?;
                self.equilibrate(at, &mut y)?;
                self.push_history(at, &y);
                self.int.restart(&self.sys, at, &y);
            } else if self.next_exec() <= t1 + 1e-12 * (1.0 + t1.abs()) {
                let mut y = y1;
                if self.handle_events(t1, &mut y)? {
                    self.equilibrate(t1, &mut y)?;
                    self.push_history(t1, &y);
                    self.int.restart(&self.sys, t1, &y);
                }
            }
        }
        Ok(())
    }

    fn push_history(&self, t: f64, y: &[f64]) {
        if self.model.has_delay {
            let mut f = vec![0.0; y.len()];
            self.sys.rhs(t, y, &mut f);
            let vals = self.sys.values(t, y);
            self.sys.history.borrow_mut().push(HistoryPoint { t, y: y.to_vec(), f, vals });
        }
    }
}

/// The fast reactions alone, in pseudo-time at a frozen simulation time.
struct FastOnly<'s, 'a> {
    sys: &'s KineticSystem<'a>,
    t: f64,
    mass: Vec<f64>,
}

impl OdeSystem for FastOnly<'_, '_> {
    fn dim(&self) -> usize {
        self.sys.dim()
    }

    fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        let vals = self.sys.fill_vals(self.t, y);
        self.sys.model.fast_rhs_from_vals(&vals, dy);
        self.sys.put_vals(vals);
    }

    fn jacobian(&self, t: f64, y: &[f64], jac: &mut DMatrix<f64>) {
        let n = y.len();
        let mut base = vec![0.0; n];
        self.rhs(t, y, &mut base);
        let mut plus = vec![0.0; n];
        let mut probe = y.to_vec();
        for j in 0..n {
            let h = (f64::EPSILON * y[j].abs().max(1e-5)).sqrt();
            probe[j] = y[j] + h;
            self.rhs(t, &probe, &mut plus);
            probe[j] = y[j];
            for i in 0..n {
                jac[(i, j)] = (plus[i] - base[i]) / h;
            }
        }
    }

    fn time_derivative(&self, _t: f64, _y: &[f64], out: &mut [f64]) {
        out.fill(0.0);
    }

    fn mass(&self) -> &[f64] {
        &self.mass
    }
}

/// Simulates `model` and records all slots at each of `times`, which must not decrease.
pub fn simulate(model: &KineticModel, times: &[f64], options: &SimOptions) -> Result<TimeCourse> {
    if times.is_empty() {
        return Err(Error::Invalid("at least one output time is required".into()));
    }
    if times.windows(2).any(|w| w[1] < w[0]) {
        return Err(Error::Invalid("output times must not decrease".into()));
    }
    if times[0] < 0.0 {
        return Err(Error::Invalid("output times must not be negative".into()));
    }
    // The dynamics start at time zero whatever the first output time.
    let t0 = 0.0;
    let mut y = model.initial_state();
    let sys = KineticSystem::new(model, options.numeric_jacobian);
    let int: Box<dyn Integrator> = match options.method {
        Method::Rodas4 => Box::new(Rodas4::new()),
        Method::Bdf => Box::new(Bdf::new()),
    };
    let mut h_max = f64::INFINITY;
    if model.has_delay {
        h_max = model.min_constant_lag.unwrap_or(0.1).min(0.5);
    }
    let tol = Tolerances { rtol: options.rtol, atol: options.atol, h_max, ..Tolerances::default() };
    let mut driver = Driver {
        model,
        sys,
        int,
        tol,
        trig: model.events.iter().map(|e| e.initial_value).collect(),
        pending: Vec::new(),
        rng: ChaCha8Rng::seed_from_u64(options.seed),
        rate_thresholds: Vec::new(),
        trigger_thresholds: Vec::new(),
        steps: model.time_steps(),
        end: times[times.len() - 1],
        event_times: Vec::new(),
    };
    driver.sys.t0 = t0;
    driver.find_thresholds(t0, &y)?;
    driver.equilibrate(t0, &mut y)?;
    driver.handle_events(t0, &mut y)?;
    driver.equilibrate(t0, &mut y)?;
    driver.sys.y0.clone_from(&y);
    driver.sys.vals0 = driver.sys.values(t0, &y);
    driver.push_history(t0, &y);
    driver.int.restart(&driver.sys, t0, &y);
    let mut out = TimeCourse {
        time: Vec::with_capacity(times.len()),
        values: Vec::with_capacity(times.len()),
        stats: Stats::default(),
        event_times: Vec::new(),
        pending_until: None,
    };
    for &t in times {
        driver.advance(t)?;
        out.time.push(t);
        out.values.push(driver.sys.values(t, driver.int.state()));
    }
    out.stats = driver.int.stats();
    out.pending_until =
        driver.pending.iter().map(|p| p.exec_time).fold(None, |a: Option<f64>, b| Some(a.map_or(b, |a| a.max(b))));
    out.event_times = driver.event_times;
    Ok(out)
}
