//! Stochastic simulation of reaction networks at low copy numbers.
//!
//! [`Algorithm::Direct`] is Gillespie's exact direct method (1977). [`Algorithm::TauLeap`] is the
//! adaptive explicit tau-leaping of Cao, Gillespie and Petzold (2006), with Poisson increments and a
//! fallback to direct steps when few reactions are critical; it is an approximation whose error is
//! set by `epsilon`. Rate rules on non-species symbols and events integrate with the exact trajectory
//! between firings by a piecewise-constant approximation of the step, so models that use them are
//! exact only as the propensities change slowly relative to the time between firings.

use super::compile::{Coef, KineticModel, NoDelay, eval};
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

struct Network {
    /// Per reaction: `(state, change)`.
    changes: Vec<Vec<(usize, f64)>>,
}

fn network(model: &KineticModel) -> Result<Network> {
    let mut changes: Vec<Vec<(usize, f64)>> = vec![Vec::new(); model.reactions.len()];
    for e in &model.stoich {
        let coef = match e.coef {
            Coef::Const(v) => v,
            Coef::Slot(s) => model.initial[s],
        };
        let conv = e.conversion.map_or(1.0, |s| model.initial[s]);
        match changes[e.reaction].iter_mut().find(|(s, _)| *s == e.state) {
            Some(entry) => entry.1 += e.sign * coef * conv,
            None => changes[e.reaction].push((e.state, e.sign * coef * conv)),
        }
    }
    Ok(Network { changes })
}

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

struct Trajectory<'a> {
    model: &'a KineticModel,
    net: &'a Network,
    y: Vec<f64>,
    vals: Vec<f64>,
    t: f64,
    triggers: Vec<bool>,
    pending: Vec<(usize, f64, Option<Vec<f64>>)>,
}

impl Trajectory<'_> {
    fn refresh(&mut self) {
        self.model.eval_state(self.t, &self.y, &mut self.vals, &NoDelay);
    }

    fn propensities(&self, out: &mut [f64]) -> f64 {
        let mut total = 0.0;
        for (j, r) in self.model.reactions.iter().enumerate() {
            let a = self.vals[r.rate].max(0.0);
            out[j] = if a.is_finite() { a } else { 0.0 };
            total += out[j];
        }
        total
    }

    /// Advances rate rules by `dt` at the current values.
    fn advance_rules(&mut self, dt: f64) {
        for (i, c) in &self.model.rate_rules {
            let d = eval(c, &self.vals, self.t, &NoDelay);
            self.y[*i] += d * dt;
        }
    }

    fn fire(&mut self, j: usize, count: f64) {
        for &(state, change) in &self.net.changes[j] {
            self.y[state] += change * count;
        }
    }

    fn check_events(&mut self) {
        if self.model.events.is_empty() {
            return;
        }
        for _ in 0..1000 {
            self.refresh();
            let current: Vec<bool> =
                self.model.events.iter().map(|e| eval(&e.trigger, &self.vals, self.t, &NoDelay) != 0.0).collect();
            for (k, e) in self.model.events.iter().enumerate() {
                if current[k] && !self.triggers[k] {
                    let delay = e.delay.as_ref().map_or(0.0, |d| eval(d, &self.vals, self.t, &NoDelay));
                    let values = e
                        .use_values_from_trigger_time
                        .then(|| e.assignments.iter().map(|(_, c)| eval(c, &self.vals, self.t, &NoDelay)).collect());
                    self.pending.push((k, self.t + delay, values));
                } else if !current[k] && self.triggers[k] && !e.persistent {
                    self.pending.retain(|p| p.0 != k);
                }
            }
            self.triggers = current;
            let Some(pos) = self.pending.iter().position(|p| p.1 <= self.t) else { return };
            let (k, _, values) = self.pending.remove(pos);
            let e = &self.model.events[k];
            let values: Vec<f64> = values
                .unwrap_or_else(|| e.assignments.iter().map(|(_, c)| eval(c, &self.vals, self.t, &NoDelay)).collect());
            for ((target, _), value) in e.assignments.iter().zip(values) {
                match *target {
                    super::compile::Target::State(i) => self.y[i] = value,
                    super::compile::Target::Amount { state, compartment, hosu } => {
                        self.y[state] = if hosu { value } else { value * self.vals[compartment] };
                    }
                }
            }
        }
    }
}

/// Largest step for tau-leaping (Cao, Gillespie and Petzold 2006, equation 8), treating all
/// reactant species as of first order for the highest-order bound `g = 2`.
fn leap_size(model: &KineticModel, net: &Network, y: &[f64], a: &[f64], epsilon: f64) -> f64 {
    let n = y.len();
    let mut mu = vec![0.0; n];
    let mut sigma = vec![0.0; n];
    for (j, changes) in net.changes.iter().enumerate() {
        for &(s, c) in changes {
            mu[s] += c * a[j];
            sigma[s] += c * c * a[j];
        }
    }
    let mut tau = f64::INFINITY;
    let _ = model;
    for s in 0..n {
        if sigma[s] == 0.0 {
            continue;
        }
        let bound = (epsilon * y[s].abs() / 2.0).max(1.0);
        tau = tau.min(bound / mu[s].abs().max(1e-300)).min(bound * bound / sigma[s]);
    }
    tau
}

fn realize(
    model: &KineticModel,
    net: &Network,
    times: &[f64],
    algorithm: Algorithm,
    seed: u64,
    out: &mut [Vec<f64>],
) -> Result<()> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut traj = Trajectory {
        model,
        net,
        y: model.initial_state(),
        vals: model.initial.clone(),
        t: times[0],
        triggers: model.events.iter().map(|e| e.initial_value).collect(),
        pending: Vec::new(),
    };
    let mut a = vec![0.0; model.reactions.len()];
    traj.check_events();
    let mut next_out = 0;
    let t_end = times[times.len() - 1];
    let has_rules = !model.rate_rules.is_empty();
    // Rate rules and events change propensities between firings, so steps are bounded.
    let bound = if has_rules || !model.events.is_empty() { (t_end - times[0]) / 2000.0 } else { f64::INFINITY };
    loop {
        traj.refresh();
        let total = traj.propensities(&mut a);
        let (dt, leap) = match algorithm {
            Algorithm::Direct => {
                let dt = if total > 0.0 { -(1.0 - rng.random::<f64>()).ln() / total } else { f64::INFINITY };
                (dt, false)
            }
            Algorithm::TauLeap { epsilon } => {
                let tau = leap_size(model, net, &traj.y, &a, epsilon);
                if total > 0.0 && tau * total > 10.0 {
                    (tau, true)
                } else {
                    let dt = if total > 0.0 { -(1.0 - rng.random::<f64>()).ln() / total } else { f64::INFINITY };
                    (dt, false)
                }
            }
        };
        let next_pending = traj.pending.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let stop_at = times[next_out].min(next_pending).min(traj.t + bound);
        if traj.t + dt >= stop_at {
            // No firing before the next stop.
            traj.advance_rules(stop_at - traj.t);
            traj.t = stop_at;
            traj.check_events();
            while next_out < times.len() && times[next_out] <= traj.t {
                traj.refresh();
                out[next_out].copy_from_slice(&traj.vals);
                next_out += 1;
            }
            if next_out >= times.len() {
                return Ok(());
            }
            continue;
        }
        traj.advance_rules(dt);
        traj.t += dt;
        if leap {
            // Poisson increments, halving the leap when an amount would turn negative.
            let snapshot = traj.y.clone();
            let mut tau = dt;
            loop {
                for (j, &rate) in a.iter().enumerate() {
                    let count = poisson(&mut rng, rate * tau);
                    traj.fire(j, count);
                }
                if traj.y.iter().all(|&v| v >= -1e-9) {
                    break;
                }
                traj.y.clone_from(&snapshot);
                tau *= 0.5;
            }
        } else {
            let mut pick = rng.random::<f64>() * total;
            let mut chosen = a.len() - 1;
            for (j, &rate) in a.iter().enumerate() {
                if pick < rate {
                    chosen = j;
                    break;
                }
                pick -= rate;
            }
            traj.fire(chosen, 1.0);
        }
        traj.check_events();
    }
}

/// Running mean, sum of squared deviations and count, as accumulated by one worker.
type Accumulated = (Vec<Vec<f64>>, Vec<Vec<f64>>, f64);

/// Runs `runs` realizations in parallel and returns per-slot mean and standard deviation.
pub fn ensemble(model: &KineticModel, times: &[f64], runs: usize, algorithm: Algorithm, seed: u64) -> Result<Ensemble> {
    if times.is_empty() || runs == 0 {
        return Err(Error::Invalid("an ensemble needs output times and at least one run".into()));
    }
    let net = network(model)?;
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
            realize(
                model,
                &net,
                times,
                algorithm,
                seed.wrapping_add(id.wrapping_mul(0x9E37_79B9_7F4A_7C15)),
                &mut sample,
            )?;
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
