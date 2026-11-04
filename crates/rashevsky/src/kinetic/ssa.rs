//! Stochastic simulation of reaction networks at low copy numbers.
//!
//! [`Algorithm::Direct`] is Gillespie's direct method (1977), exact while propensities stay constant
//! between firings. [`Algorithm::TauLeap`] is the explicit tau-leaping of Cao, Gillespie and Petzold
//! (2006), an approximation whose error `epsilon` sets, with critical reactions and the step selection
//! of their equations 32 and 33. Reaction orders for that selection come from the reactant
//! stoichiometry, which is the true order for mass-action laws.
//!
//! A propensity that depends on the time or on a rate rule is held constant over a step of at most
//! 1/2000 of the run, which makes such models approximate. Models outside what a propensity can
//! describe are refused: negative or undefined rates, algebraic rules, `delay` and fast reactions.
//! A reaction whose firing would leave a species below zero has propensity zero, whatever its law says.

use super::compile::{Coef, KineticModel, NoDelay, Target, eval};
use crate::error::{Error, Result};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};

/// Stochastic algorithm.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum Algorithm {
    /// Exact direct method.
    Direct,
    /// Adaptive explicit tau-leaping.
    TauLeap {
        /// Bound on the relative change of propensities per leap.
        epsilon: f64,
    },
}

/// Mean and standard deviation over an ensemble.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Ensemble {
    /// Output times.
    pub time: Vec<f64>,
    /// `mean[k][slot]`.
    pub mean: Vec<Vec<f64>>,
    /// `sd[k][slot]`.
    pub sd: Vec<Vec<f64>>,
    /// Number of realizations.
    pub runs: usize,
}

/// A reaction is critical when fewer than this many firings would use up a reactant (Cao 2006, step 1).
const CRITICAL_FIRINGS: f64 = 10.0;
/// Below this many expected firings a leap gains nothing over direct steps (Cao 2006, step 3).
const LEAP_MINIMUM: f64 = 10.0;
/// Direct steps taken before a leap is considered again (Cao 2006, step 3).
const DIRECT_BURST: usize = 100;

fn poisson(rng: &mut ChaCha8Rng, mean: f64) -> f64 {
    if mean <= 0.0 {
        return 0.0;
    }
    if mean < 30.0 {
        // Knuth's product method.
        let limit = (-mean).exp();
        let mut k = 0.0;
        let mut p = 1.0;
        loop {
            p *= rng.random::<f64>();
            if p <= limit {
                return k;
            }
            k += 1.0;
        }
    }
    // Atkinson's rejection method for large means.
    let c = 0.767 - 3.36 / mean;
    let beta = std::f64::consts::PI / (3.0 * mean).sqrt();
    let alpha = beta * mean;
    let k = c.ln() - mean - beta.ln();
    loop {
        let u: f64 = rng.random();
        if u <= 0.0 || u >= 1.0 {
            continue;
        }
        let x = (alpha - ((1.0 - u) / u).ln()) / beta;
        let n = (x + 0.5).floor();
        if n < 0.0 {
            continue;
        }
        let v: f64 = rng.random();
        if v <= 0.0 {
            continue;
        }
        let y = alpha - beta * x;
        let lhs = y + (v / (1.0 + y.exp()).powi(2)).ln();
        let rhs = k + n * mean.ln() - ln_gamma(n + 1.0);
        if lhs <= rhs {
            return n;
        }
    }
}

/// Lanczos approximation of the logarithm of the gamma function.
fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    let z = x - 1.0;
    let mut a = C[0];
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (z + i as f64);
    }
    let t = z + G + 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + a.ln()
}

/// Reactants of a reaction at the current coefficients, with the molecules one firing takes.
struct Reactant {
    state: usize,
    taken: f64,
}

struct Trajectory<'a> {
    model: &'a KineticModel,
    y: Vec<f64>,
    vals: Vec<f64>,
    t: f64,
    triggers: Vec<bool>,
    pending: Vec<(usize, f64, Option<Vec<f64>>)>,
    /// Per reaction: `(state, change per firing)`.
    changes: Vec<Vec<(usize, f64)>>,
    variable_stoichiometry: bool,
}

impl Trajectory<'_> {
    fn refresh(&mut self) {
        self.model.eval_state(self.t, &self.y, &mut self.vals, &NoDelay);
        if self.variable_stoichiometry {
            self.read_changes();
        }
    }

    fn read_changes(&mut self) {
        for list in &mut self.changes {
            list.clear();
        }
        for e in &self.model.stoich {
            let coef = match e.coef {
                Coef::Const(v) => v,
                Coef::Slot(s) => self.vals[s],
            };
            let conversion = e.conversion.map_or(1.0, |s| self.vals[s]);
            let list = &mut self.changes[e.reaction];
            match list.iter_mut().find(|(s, _)| *s == e.state) {
                Some(entry) => entry.1 += e.sign * coef * conversion,
                None => list.push((e.state, e.sign * coef * conversion)),
            }
        }
    }

    fn propensities(&self, out: &mut [f64]) -> Result<f64> {
        let mut total = 0.0;
        for (j, r) in self.model.reactions.iter().enumerate() {
            let a = self.vals[r.rate];
            if a.is_nan() || a < 0.0 || a.is_infinite() {
                return Err(Error::Kinetic(format!(
                    "reaction {} has propensity {a} at t = {}; stochastic simulation needs each direction written as its own reaction",
                    r.id, self.t
                )));
            }
            let possible = self.changes[j].iter().all(|&(s, c)| c >= 0.0 || self.y[s] + c >= -1e-9);
            out[j] = if possible { a } else { 0.0 };
            total += out[j];
        }
        Ok(total)
    }

    /// Advances rate rules by `dt` at the current values.
    fn advance_rules(&mut self, dt: f64) {
        for (i, c) in &self.model.rate_rules {
            let d = eval(c, &self.vals, self.t, &NoDelay);
            self.y[*i] += d * dt;
        }
    }

    fn fire(&mut self, j: usize, count: f64) {
        for &(state, change) in &self.changes[j] {
            self.y[state] += change * count;
        }
    }

    /// Fires every trigger transition and executes what is due now, highest priority first and ties at
    /// random, re-evaluating the triggers after each execution as SBML specifies.
    fn check_events(&mut self, rng: &mut ChaCha8Rng) -> Result<()> {
        if self.model.events.is_empty() {
            return Ok(());
        }
        for _ in 0..1000 {
            self.refresh();
            let current: Vec<bool> = self
                .model
                .events
                .iter()
                .map(|e| {
                    let v = eval(&e.trigger, &self.vals, self.t, &NoDelay);
                    v != 0.0 && !v.is_nan()
                })
                .collect();
            for (k, e) in self.model.events.iter().enumerate() {
                if current[k] && !self.triggers[k] {
                    let delay = e.delay.as_ref().map_or(0.0, |d| eval(d, &self.vals, self.t, &NoDelay));
                    if delay < 0.0 || delay.is_nan() {
                        return Err(Error::Kinetic(format!("event {} has delay {delay} at t = {}", e.id, self.t)));
                    }
                    let values = e
                        .use_values_from_trigger_time
                        .then(|| e.assignments.iter().map(|(_, c)| eval(c, &self.vals, self.t, &NoDelay)).collect());
                    self.pending.push((k, self.t + delay, values));
                } else if !current[k] && self.triggers[k] && !e.persistent {
                    self.pending.retain(|p| p.0 != k);
                }
            }
            self.triggers = current;
            let mut best: Option<(usize, f64)> = None;
            let mut ties = 1u32;
            for (i, p) in self.pending.iter().enumerate() {
                if p.1 <= self.t {
                    let e = &self.model.events[p.0];
                    let priority = e.priority.as_ref().map_or(f64::NEG_INFINITY, |c| eval(c, &self.vals, self.t, &NoDelay));
                    match best {
                        None => best = Some((i, priority)),
                        Some((_, b)) if priority > b => {
                            best = Some((i, priority));
                            ties = 1;
                        }
                        Some((_, b)) if priority == b => {
                            ties += 1;
                            if rng.random_range(0..ties) == 0 {
                                best = Some((i, priority));
                            }
                        }
                        _ => {}
                    }
                }
            }
            let Some((pos, _)) = best else { return Ok(()) };
            let (k, _, values) = self.pending.remove(pos);
            let e = &self.model.events[k];
            let values: Vec<f64> = values
                .unwrap_or_else(|| e.assignments.iter().map(|(_, c)| eval(c, &self.vals, self.t, &NoDelay)).collect());
            for ((target, _), value) in e.assignments.iter().zip(values) {
                match *target {
                    Target::State(i) => self.y[i] = value,
                    Target::Amount { state, compartment, hosu } => {
                        self.y[state] = if hosu { value } else { value * self.vals[compartment] };
                    }
                }
            }
        }
        Err(Error::Kinetic(format!("the events did not settle at t = {}", self.t)))
    }

    /// Reactants of each reaction, from the current change vectors.
    fn reactants(&self) -> Vec<Vec<Reactant>> {
        self.changes
            .iter()
            .map(|list| list.iter().filter(|(_, c)| *c < 0.0).map(|&(state, c)| Reactant { state, taken: -c }).collect())
            .collect()
    }
}

/// `g_i` of Cao et al. (2006), equation 27: the bound on a propensity's relative change per unit
/// relative change of species `i`, set by the highest-order reaction taking `i` as a reactant.
fn relative_change_bound(order: f64, taken: f64, x: f64) -> f64 {
    let lone = |v: f64| if v > 1.0 { 1.0 / v } else { 0.0 };
    if order <= 1.0 {
        1.0
    } else if order <= 2.0 {
        if taken >= 2.0 { 2.0 + lone(x - 1.0) } else { 2.0 }
    } else if taken >= 3.0 {
        3.0 + lone(x - 1.0) + 2.0 * lone(x - 2.0)
    } else if taken >= 2.0 {
        1.5 * (2.0 + lone(x - 1.0))
    } else {
        3.0
    }
}

/// The leap `tau'` of equation 33, over non-critical reactions and their reactant species.
fn leap_size(y: &[f64], reactants: &[Vec<Reactant>], changes: &[Vec<(usize, f64)>], a: &[f64], critical: &[bool], epsilon: f64) -> f64 {
    let n = y.len();
    let mut mu = vec![0.0; n];
    let mut sigma = vec![0.0; n];
    for (j, list) in changes.iter().enumerate() {
        if critical[j] || a[j] == 0.0 {
            continue;
        }
        for &(s, c) in list {
            mu[s] += c * a[j];
            sigma[s] += c * c * a[j];
        }
    }
    let mut highest = vec![(0.0f64, 0.0f64); n];
    for list in reactants {
        let order: f64 = list.iter().map(|r| r.taken).sum();
        for r in list {
            let entry = &mut highest[r.state];
            if order > entry.0 || (order == entry.0 && r.taken > entry.1) {
                *entry = (order, r.taken);
            }
        }
    }
    let mut tau = f64::INFINITY;
    for s in 0..n {
        let (order, taken) = highest[s];
        if order == 0.0 || sigma[s] == 0.0 {
            continue;
        }
        let g = relative_change_bound(order, taken, y[s]);
        let bound = (epsilon * y[s].abs() / g).max(1.0);
        tau = tau.min(bound / mu[s].abs().max(1e-300)).min(bound * bound / sigma[s]);
    }
    tau
}

/// Settings of one realization.
struct Plan<'a> {
    times: &'a [f64],
    algorithm: Algorithm,
    thresholds: &'a [f64],
    /// Longest step over which propensities are held constant.
    bound: f64,
}

/// Records every output time at or before the current time; true once all of them are recorded.
fn record(traj: &mut Trajectory<'_>, times: &[f64], next: &mut usize, out: &mut [Vec<f64>]) -> bool {
    while *next < times.len() && times[*next] <= traj.t {
        traj.refresh();
        out[*next].copy_from_slice(&traj.vals);
        *next += 1;
    }
    *next >= times.len()
}

fn direct_step(traj: &mut Trajectory<'_>, a: &[f64], total: f64, stop_at: f64, rng: &mut ChaCha8Rng) {
    let dt = if total > 0.0 { -(1.0 - rng.random::<f64>()).ln() / total } else { f64::INFINITY };
    if traj.t + dt >= stop_at {
        // The waiting time is memoryless, so nothing fires before the stop.
        traj.advance_rules(stop_at - traj.t);
        traj.t = stop_at;
        return;
    }
    traj.advance_rules(dt);
    traj.t += dt;
    let mut pick = rng.random::<f64>() * total;
    let mut chosen = None;
    for (j, &rate) in a.iter().enumerate() {
        if rate > 0.0 {
            chosen = Some(j);
            if pick < rate {
                break;
            }
            pick -= rate;
        }
    }
    if let Some(j) = chosen {
        traj.fire(j, 1.0);
    }
}

/// One leap by steps 2 to 6 of Cao et al. (2006). Returns false when direct steps should be used.
fn leap(traj: &mut Trajectory<'_>, a: &[f64], total: f64, stop_at: f64, epsilon: f64, rng: &mut ChaCha8Rng) -> bool {
    let reactants = traj.reactants();
    let critical: Vec<bool> = reactants
        .iter()
        .zip(a)
        .map(|(list, &rate)| {
            rate > 0.0 && list.iter().map(|r| (traj.y[r.state] / r.taken).floor()).fold(f64::INFINITY, f64::min) < CRITICAL_FIRINGS
        })
        .collect();
    let mut tau_noncritical = leap_size(&traj.y, &reactants, &traj.changes, a, &critical, epsilon);
    let critical_total: f64 = a.iter().zip(&critical).filter(|&(_, c)| *c).map(|(r, _)| r).sum();
    let snapshot = traj.y.clone();
    loop {
        if tau_noncritical * total < LEAP_MINIMUM {
            return false;
        }
        let tau_critical = if critical_total > 0.0 { -(1.0 - rng.random::<f64>()).ln() / critical_total } else { f64::INFINITY };
        let remaining = stop_at - traj.t;
        let tau = tau_noncritical.min(tau_critical).min(remaining);
        for (j, &rate) in a.iter().enumerate() {
            if !critical[j] && rate > 0.0 {
                let count = poisson(rng, rate * tau);
                traj.fire(j, count);
            }
        }
        if tau_critical <= tau {
            let mut pick = rng.random::<f64>() * critical_total;
            let mut chosen = None;
            for (j, &rate) in a.iter().enumerate() {
                if critical[j] && rate > 0.0 {
                    chosen = Some(j);
                    if pick < rate {
                        break;
                    }
                    pick -= rate;
                }
            }
            if let Some(j) = chosen {
                traj.fire(j, 1.0);
            }
        }
        // Only species are bound below by zero; parameters with rate rules may be negative.
        let negative = reactants.iter().flatten().any(|r| traj.y[r.state] < -1e-9);
        if !negative {
            traj.advance_rules(tau);
            traj.t += tau;
            return true;
        }
        traj.y.clone_from(&snapshot);
        tau_noncritical *= 0.5;
    }
}

fn realize(model: &KineticModel, plan: &Plan<'_>, seed: u64, out: &mut [Vec<f64>]) -> Result<()> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut traj = Trajectory {
        model,
        y: model.initial_state(),
        vals: model.initial.clone(),
        t: 0.0,
        triggers: model.events.iter().map(|e| e.initial_value).collect(),
        pending: Vec::new(),
        changes: vec![Vec::new(); model.reactions.len()],
        variable_stoichiometry: model.has_varying_stoichiometry(),
    };
    traj.vals.clone_from(&model.initial);
    traj.read_changes();
    let mut a = vec![0.0; model.reactions.len()];
    traj.check_events(&mut rng)?;
    let mut next_out = 0;
    if record(&mut traj, plan.times, &mut next_out, out) {
        return Ok(());
    }
    let mut direct_left = 0usize;
    loop {
        traj.refresh();
        let total = traj.propensities(&mut a)?;
        let next_pending = traj.pending.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let threshold = plan.thresholds.iter().copied().find(|&t| t > traj.t).unwrap_or(f64::INFINITY);
        let stop_at = plan.times[next_out].min(next_pending).min(threshold).min(traj.t + plan.bound);
        match plan.algorithm {
            Algorithm::Direct => direct_step(&mut traj, &a, total, stop_at, &mut rng),
            Algorithm::TauLeap { epsilon } => {
                if direct_left > 0 {
                    direct_left -= 1;
                    direct_step(&mut traj, &a, total, stop_at, &mut rng);
                } else if total == 0.0 || !leap(&mut traj, &a, total, stop_at, epsilon, &mut rng) {
                    direct_left = DIRECT_BURST - 1;
                    direct_step(&mut traj, &a, total, stop_at, &mut rng);
                }
            }
        }
        traj.check_events(&mut rng)?;
        if record(&mut traj, plan.times, &mut next_out, out) {
            return Ok(());
        }
    }
}

/// Running mean, sum of squared deviations and count, as accumulated by one worker.
type Accumulated = (Vec<Vec<f64>>, Vec<Vec<f64>>, f64);

/// Runs `runs` realizations in parallel and returns per-slot mean and standard deviation.
pub fn ensemble(model: &KineticModel, times: &[f64], runs: usize, algorithm: Algorithm, seed: u64) -> Result<Ensemble> {
    if times.is_empty() || runs == 0 {
        return Err(Error::Invalid("an ensemble needs output times and at least one run".into()));
    }
    if times[0] < 0.0 || times.windows(2).any(|w| w[1] < w[0]) {
        return Err(Error::Invalid("output times must not be negative or decrease".into()));
    }
    if !model.algebraic.is_empty() {
        return Err(Error::Kinetic("stochastic simulation does not cover algebraic rules".into()));
    }
    if model.has_delay {
        return Err(Error::Kinetic("stochastic simulation does not cover delay".into()));
    }
    if model.has_fast() {
        return Err(Error::Kinetic("stochastic simulation does not cover fast reactions".into()));
    }
    if let Algorithm::TauLeap { epsilon } = algorithm
        && !(epsilon > 0.0 && epsilon < 1.0)
    {
        return Err(Error::Invalid("the tau-leaping epsilon must lie between 0 and 1".into()));
    }
    let end = times[times.len() - 1];
    let continuous = model.reads_time() || !model.rate_rules.is_empty();
    let mut thresholds = model.rate_thresholds();
    thresholds.extend(model.trigger_thresholds());
    thresholds.sort_by(f64::total_cmp);
    thresholds.dedup();
    let plan = Plan {
        times,
        algorithm,
        thresholds: &thresholds,
        bound: if continuous && end > 0.0 { end / 2000.0 } else { f64::INFINITY },
    };
    let slots = model.symbols.len();
    let ids: Vec<u64> = (0..runs as u64).collect();
    let chunk = crate::parallel::chunk_size(runs).max(1);
    // Welford accumulation per chunk, merged at the end.
    let parts: Vec<Result<Accumulated>> = crate::parallel::map_chunks(&ids, chunk, |slice| {
        let mut mean = vec![vec![0.0; slots]; times.len()];
        let mut m2 = vec![vec![0.0; slots]; times.len()];
        let mut sample = vec![vec![0.0; slots]; times.len()];
        let mut count = 0.0;
        for &id in slice {
            realize(model, &plan, seed.wrapping_add(id.wrapping_mul(0x9E37_79B9_7F4A_7C15)), &mut sample)?;
            count += 1.0;
            for k in 0..times.len() {
                for s in 0..slots {
                    let delta = sample[k][s] - mean[k][s];
                    mean[k][s] += delta / count;
                    m2[k][s] += delta * (sample[k][s] - mean[k][s]);
                }
            }
        }
        Ok((mean, m2, count))
    });
    let mut mean = vec![vec![0.0; slots]; times.len()];
    let mut m2 = vec![vec![0.0; slots]; times.len()];
    let mut total = 0.0;
    for part in parts {
        let (pm, pm2, pc) = part?;
        for k in 0..times.len() {
            for s in 0..slots {
                let delta = pm[k][s] - mean[k][s];
                let new_total = total + pc;
                mean[k][s] += delta * pc / new_total;
                m2[k][s] += pm2[k][s] + delta * delta * total * pc / new_total;
            }
        }
        total += pc;
    }
    let sd = m2.iter().map(|row| row.iter().map(|v| (v / (total - 1.0).max(1.0)).sqrt()).collect()).collect();
    Ok(Ensemble { time: times.to_vec(), mean, sd, runs })
}
