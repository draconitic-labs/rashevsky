//! Stiff integrators for `M y' = f(t, y)` with a diagonal mass matrix.
//!
//! [`Rodas4`] is the L-stable, stiffly accurate Rosenbrock method of order 4 with an embedded order-3
//! error estimate (Hairer and Wanner, *Solving Ordinary Differential Equations II*, coefficients from
//! their `rodas.f`). [`Bdf`] is a variable-order (1 to 5), quasi-constant-step backward difference
//! method in the numerical differentiation form of Shampine and Reichelt (1997). Entries of zero in
//! the mass diagonal mark algebraic equations, both methods handle index-1 systems.

use crate::error::{Error, Result};
use nalgebra::{DMatrix, DVector};

/// The right-hand side and its derivatives.
pub trait OdeSystem {
    /// Number of equations.
    fn dim(&self) -> usize;
    /// Evaluates `f(t, y)` into `dy`.
    fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]);
    /// Evaluates the Jacobian `df/dy` into `jac`, which is zero on entry.
    fn jacobian(&self, t: f64, y: &[f64], jac: &mut DMatrix<f64>);
    /// Evaluates `df/dt` into `out`, which is zero on entry.
    fn time_derivative(&self, t: f64, y: &[f64], out: &mut [f64]);
    /// Diagonal of the mass matrix: 1 for differential equations, 0 for algebraic ones.
    fn mass(&self) -> &[f64];
}

/// Tolerances and limits of an integration.
#[derive(Clone, Copy, Debug)]
pub struct Tolerances {
    /// Relative tolerance.
    pub rtol: f64,
    /// Absolute tolerance.
    pub atol: f64,
    /// Largest step allowed.
    pub h_max: f64,
    /// Largest number of attempted steps per call to `advance`.
    pub max_steps: usize,
}

impl Default for Tolerances {
    fn default() -> Self {
        Self { rtol: 1e-8, atol: 1e-12, h_max: f64::INFINITY, max_steps: 2_000_000 }
    }
}

/// Counters of work done by an integrator.
#[derive(Clone, Copy, Debug, Default)]
pub struct Stats {
    /// Accepted steps.
    pub steps: u64,
    /// Rejected steps.
    pub rejected: u64,
    /// Right-hand side evaluations.
    pub rhs_evals: u64,
    /// Jacobian evaluations.
    pub jacobians: u64,
    /// Matrix factorizations.
    pub factorizations: u64,
}

/// A one-step-at-a-time integrator that can stop at a given time and be restarted after a jump.
pub trait Integrator {
    /// Restarts from `(t, y)`, discarding step history.
    fn restart(&mut self, sys: &dyn OdeSystem, t: f64, y: &[f64]);
    /// Current time.
    fn time(&self) -> f64;
    /// Current state.
    fn state(&self) -> &[f64];
    /// Takes one accepted step that does not pass `t_limit`.
    fn step(&mut self, sys: &dyn OdeSystem, t_limit: f64, tol: &Tolerances) -> Result<()>;
    /// Interpolates inside the last accepted step into `out`.
    fn interpolate(&self, t: f64, out: &mut [f64]);
    /// Start of the last accepted step.
    fn previous_time(&self) -> f64;
    /// Work counters.
    fn stats(&self) -> Stats;
}

fn wrms(values: &[f64], y0: &[f64], y1: &[f64], tol: &Tolerances) -> f64 {
    let n = values.len().max(1) as f64;
    let mut sum = 0.0;
    for k in 0..values.len() {
        let sk = tol.atol + tol.rtol * y0[k].abs().max(y1[k].abs());
        let r = values[k] / sk;
        sum += r * r;
    }
    (sum / n).sqrt()
}

fn all_finite(v: &[f64]) -> bool {
    v.iter().all(|x| x.is_finite())
}

/// Initial step following Hairer, Norsett and Wanner, section II.4.
fn initial_step(sys: &dyn OdeSystem, t: f64, y: &[f64], f0: &[f64], tol: &Tolerances, order: f64) -> f64 {
    let n = y.len();
    if n == 0 {
        return 1e-3;
    }
    let mass = sys.mass();
    let scale: Vec<f64> = y.iter().map(|v| tol.atol + tol.rtol * v.abs()).collect();
    let norm =
        |v: &[f64]| -> f64 { (v.iter().zip(&scale).map(|(a, s)| (a / s) * (a / s)).sum::<f64>() / n as f64).sqrt() };
    let d0 = norm(y);
    let diff: Vec<f64> = f0.iter().zip(mass).map(|(f, m)| f * m).collect();
    let d1 = norm(&diff);
    let h0 = if d0 < 1e-5 || d1 < 1e-5 { 1e-6 } else { 0.01 * d0 / d1 };
    let y1: Vec<f64> = y.iter().zip(&diff).map(|(a, f)| a + h0 * f).collect();
    let mut f1 = vec![0.0; n];
    sys.rhs(t + h0, &y1, &mut f1);
    let d2: Vec<f64> = f1.iter().zip(f0).zip(mass).map(|((a, b), m)| (a - b) * m).collect();
    let d2 = norm(&d2) / h0;
    let h1 =
        if d1 <= 1e-15 && d2 <= 1e-15 { (h0 * 1e-3).max(1e-6) } else { (0.01 / d1.max(d2)).powf(1.0 / (order + 1.0)) };
    (100.0 * h0).min(h1).min(tol.h_max)
}

// Coefficients as printed in Hairer's rodas.f, with all their digits.
#[allow(clippy::excessive_precision)]
mod coefficients {
    pub const GAMMA: f64 = 0.25;
    pub const C2: f64 = 0.386;
    pub const C3: f64 = 0.21;
    pub const C4: f64 = 0.63;
    pub const D1: f64 = 0.25;
    pub const D2: f64 = -0.1043;
    pub const D3: f64 = 0.1035;
    pub const D4: f64 = -0.3620000000000023e-1;
    pub const A21: f64 = 0.1544000000000000e1;
    pub const A31: f64 = 0.9466785280815826;
    pub const A32: f64 = 0.2557011698983284;
    pub const A41: f64 = 0.3314825187068521e1;
    pub const A42: f64 = 0.2896124015972201e1;
    pub const A43: f64 = 0.9986419139977817;
    pub const A51: f64 = 0.1221224509226641e1;
    pub const A52: f64 = 0.6019134481288629e1;
    pub const A53: f64 = 0.1253708332932087e2;
    pub const A54: f64 = -0.6878860361058950;
    pub const C21: f64 = -0.5668800000000000e1;
    pub const C31: f64 = -0.2430093356833875e1;
    pub const C32: f64 = -0.2063599157091915;
    pub const C41: f64 = -0.1073529058151375;
    pub const C42: f64 = -0.9594562251023355e1;
    pub const C43: f64 = -0.2047028614809616e2;
    pub const C51: f64 = 0.7496443313967647e1;
    pub const C52: f64 = -0.1024680431464352e2;
    pub const C53: f64 = -0.3399990352819905e2;
    pub const C54: f64 = 0.1170890893206160e2;
    pub const C61: f64 = 0.8083246795921522e1;
    pub const C62: f64 = -0.7981132988064893e1;
    pub const C63: f64 = -0.3152159432874371e2;
    pub const C64: f64 = 0.1631930543123136e2;
    pub const C65: f64 = -0.6058818238834054e1;
}
use coefficients::*;

/// Rosenbrock method RODAS of order 4(3).
#[derive(Clone, Debug, Default)]
pub struct Rodas4 {
    t: f64,
    y: Vec<f64>,
    f: Vec<f64>,
    t_prev: f64,
    y_prev: Vec<f64>,
    f_prev: Vec<f64>,
    mass: Vec<f64>,
    h: f64,
    stats: Stats,
}

impl Rodas4 {
    /// A new integrator; call [`Integrator::restart`] before stepping.
    pub fn new() -> Self {
        Self::default()
    }

    /// Takes one attempt with step `h`, returning the new state and the error norm.
    fn attempt(
        &mut self,
        sys: &dyn OdeSystem,
        h: f64,
        jac: &DMatrix<f64>,
        fx: &[f64],
        tol: &Tolerances,
    ) -> Result<Option<(Vec<f64>, f64)>> {
        let n = self.y.len();
        let mass = sys.mass();
        let fac = 1.0 / (h * GAMMA);
        let mut e = -jac.clone();
        for k in 0..n {
            e[(k, k)] += fac * mass[k];
        }
        self.stats.factorizations += 1;
        // Models driven only by events have an empty or diagonal Jacobian, which needs no factorization.
        let is_diagonal = (0..n).all(|j| (0..n).all(|i| i == j || e[(i, j)] == 0.0));
        let diagonal: Vec<f64> = (0..n).map(|k| e[(k, k)]).collect();
        let lu = if is_diagonal { None } else { Some(e.lu()) };
        let solve = |mut rhs: Vec<f64>| -> Option<Vec<f64>> {
            match &lu {
                Some(lu) => lu.solve(&DVector::from_vec(rhs)).map(|v| v.as_slice().to_vec()),
                None => {
                    for (r, d) in rhs.iter_mut().zip(&diagonal) {
                        if *d == 0.0 {
                            return None;
                        }
                        *r /= d;
                    }
                    Some(rhs)
                }
            }
        };
        let y = &self.y;
        let t = self.t;
        let mut dy = vec![0.0; n];
        let mut ynew = vec![0.0; n];
        // Stage 1.
        let rhs: Vec<f64> = (0..n).map(|k| self.f[k] + h * D1 * fx[k]).collect();
        let Some(k1) = solve(rhs) else { return Ok(None) };
        // Stage 2.
        for i in 0..n {
            ynew[i] = y[i] + A21 * k1[i];
        }
        sys.rhs(t + C2 * h, &ynew, &mut dy);
        let rhs: Vec<f64> = (0..n).map(|k| dy[k] + h * D2 * fx[k] + mass[k] * (C21 / h * k1[k])).collect();
        let Some(k2) = solve(rhs) else { return Ok(None) };
        // Stage 3.
        for i in 0..n {
            ynew[i] = y[i] + A31 * k1[i] + A32 * k2[i];
        }
        sys.rhs(t + C3 * h, &ynew, &mut dy);
        let rhs: Vec<f64> =
            (0..n).map(|k| dy[k] + h * D3 * fx[k] + mass[k] * ((C31 * k1[k] + C32 * k2[k]) / h)).collect();
        let Some(k3) = solve(rhs) else { return Ok(None) };
        // Stage 4.
        for i in 0..n {
            ynew[i] = y[i] + A41 * k1[i] + A42 * k2[i] + A43 * k3[i];
        }
        sys.rhs(t + C4 * h, &ynew, &mut dy);
        let rhs: Vec<f64> = (0..n)
            .map(|k| dy[k] + h * D4 * fx[k] + mass[k] * ((C41 * k1[k] + C42 * k2[k] + C43 * k3[k]) / h))
            .collect();
        let Some(k4) = solve(rhs) else { return Ok(None) };
        // Stage 5.
        for i in 0..n {
            ynew[i] = y[i] + A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i];
        }
        sys.rhs(t + h, &ynew, &mut dy);
        let rhs: Vec<f64> =
            (0..n).map(|k| dy[k] + mass[k] * ((C51 * k1[k] + C52 * k2[k] + C53 * k3[k] + C54 * k4[k]) / h)).collect();
        let Some(k5) = solve(rhs) else { return Ok(None) };
        // Embedded solution, then the last stage, which is also the error estimate.
        for i in 0..n {
            ynew[i] += k5[i];
        }
        sys.rhs(t + h, &ynew, &mut dy);
        let rhs: Vec<f64> = (0..n)
            .map(|k| dy[k] + mass[k] * ((C61 * k1[k] + C62 * k2[k] + C63 * k3[k] + C64 * k4[k] + C65 * k5[k]) / h))
            .collect();
        let Some(k6) = solve(rhs) else { return Ok(None) };
        for i in 0..n {
            ynew[i] += k6[i];
        }
        self.stats.rhs_evals += 5;
        if !all_finite(&ynew) {
            return Ok(None);
        }
        let err = wrms(&k6, y, &ynew, tol);
        Ok(Some((ynew, err)))
    }
}

impl Integrator for Rodas4 {
    fn restart(&mut self, sys: &dyn OdeSystem, t: f64, y: &[f64]) {
        self.t = t;
        self.t_prev = t;
        self.y = y.to_vec();
        self.y_prev = y.to_vec();
        self.f = vec![0.0; y.len()];
        sys.rhs(t, y, &mut self.f);
        self.f_prev = self.f.clone();
        self.mass = sys.mass().to_vec();
        self.stats.rhs_evals += 1;
        self.h = 0.0;
    }

    fn time(&self) -> f64 {
        self.t
    }

    fn state(&self) -> &[f64] {
        &self.y
    }

    fn previous_time(&self) -> f64 {
        self.t_prev
    }

    fn stats(&self) -> Stats {
        self.stats
    }

    fn step(&mut self, sys: &dyn OdeSystem, t_limit: f64, tol: &Tolerances) -> Result<()> {
        let n = self.y.len();
        if n == 0 {
            self.t_prev = self.t;
            self.t = t_limit;
            return Ok(());
        }
        if self.h <= 0.0 {
            self.h = initial_step(sys, self.t, &self.y, &self.f, tol, 4.0).max(1e-14);
        }
        let mut jac = DMatrix::zeros(n, n);
        sys.jacobian(self.t, &self.y, &mut jac);
        self.stats.jacobians += 1;
        let mut fx = vec![0.0; n];
        sys.time_derivative(self.t, &self.y, &mut fx);
        let remaining = t_limit - self.t;
        // A remainder of the order of the rounding in `t` is not a step.
        let rounding = 8.0 * f64::EPSILON * t_limit.abs();
        if remaining <= rounding {
            self.t_prev = self.t;
            self.t = t_limit;
            return Ok(());
        }
        let mut last_planned = None;
        for _ in 0..tol.max_steps {
            let mut h = self.h.min(tol.h_max);
            let clipped = h >= remaining * (1.0 - 1e-12) || remaining - h <= rounding;
            if clipped {
                last_planned = Some(self.h);
                h = remaining;
            }
            if h.abs() < 1e-14 * self.t.abs().max(1e-300) {
                return Err(Error::Integrator(format!("step size underflow at t = {}", self.t)));
            }
            match self.attempt(sys, h, &jac, &fx, tol)? {
                Some((ynew, err)) if err <= 1.0 => {
                    let t_new = if clipped { t_limit } else { self.t + h };
                    let mut fnew = vec![0.0; n];
                    sys.rhs(t_new, &ynew, &mut fnew);
                    self.stats.rhs_evals += 1;
                    self.stats.steps += 1;
                    self.t_prev = self.t;
                    self.y_prev = std::mem::replace(&mut self.y, ynew);
                    self.f_prev = std::mem::replace(&mut self.f, fnew);
                    self.t = t_new;
                    let factor = (0.9 * err.max(1e-10).powf(-0.25)).clamp(1.0 / 6.0, 5.0);
                    self.h = h * factor;
                    if clipped && let Some(planned) = last_planned {
                        // Keep the planned step when the limit shortened this one.
                        self.h = self.h.max(planned);
                    }
                    return Ok(());
                }
                Some((_, err)) => {
                    self.stats.rejected += 1;
                    let factor = (0.9 * err.powf(-0.25)).clamp(1.0 / 6.0, 1.0);
                    self.h = h * factor;
                }
                None => {
                    self.stats.rejected += 1;
                    self.h = h * 0.25;
                }
            }
        }
        Err(Error::Integrator(format!("too many step rejections at t = {}", self.t)))
    }

    fn interpolate(&self, t: f64, out: &mut [f64]) {
        // Cubic Hermite interpolant over the accepted step.
        let h = self.t - self.t_prev;
        if h == 0.0 {
            out.copy_from_slice(&self.y);
            return;
        }
        let s = (t - self.t_prev) / h;
        let (h00, h10) = ((1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s), s * (1.0 - s) * (1.0 - s));
        let (h01, h11) = (s * s * (3.0 - 2.0 * s), s * s * (s - 1.0));
        for k in 0..out.len() {
            out[k] = if self.mass[k] == 0.0 {
                // The right-hand side of an algebraic equation is a residual, not a slope.
                self.y_prev[k] + s * (self.y[k] - self.y_prev[k])
            } else {
                h00 * self.y_prev[k] + h10 * h * self.f_prev[k] + h01 * self.y[k] + h11 * h * self.f[k]
            };
        }
    }
}

const MAX_ORDER: usize = 5;
const NEWTON_MAXITER: usize = 4;
const KAPPA: [f64; 6] = [0.0, -0.1850, -1.0 / 9.0, -0.0823, -0.0415, 0.0];

/// Variable-order backward differences (numerical differentiation formulas).
#[derive(Clone, Debug, Default)]
pub struct Bdf {
    t: f64,
    t_prev: f64,
    order: usize,
    h: f64,
    n_equal: usize,
    d: Vec<Vec<f64>>,
    y: Vec<f64>,
    gamma: [f64; 7],
    alpha: [f64; 6],
    err_const: [f64; 7],
    jac: Option<DMatrix<f64>>,
    jac_current: bool,
    stats: Stats,
}

fn compute_r(order: usize, factor: f64) -> DMatrix<f64> {
    let n = order + 1;
    let mut m = DMatrix::zeros(n, n);
    for j in 0..n {
        m[(0, j)] = 1.0;
    }
    for i in 1..n {
        for j in 1..n {
            m[(i, j)] = (i as f64 - 1.0 - factor * j as f64) / i as f64;
        }
    }
    // Cumulative product along columns.
    for j in 0..n {
        for i in 1..n {
            m[(i, j)] *= m[(i - 1, j)];
        }
    }
    m
}

impl Bdf {
    /// A new integrator; call [`Integrator::restart`] before stepping.
    pub fn new() -> Self {
        let mut gamma = [0.0; 7];
        for k in 1..7 {
            gamma[k] = gamma[k - 1] + 1.0 / k as f64;
        }
        let mut alpha = [0.0; 6];
        for k in 0..6 {
            alpha[k] = (1.0 - KAPPA[k]) * gamma[k];
        }
        let mut err_const = [0.0; 7];
        for k in 0..6 {
            err_const[k] = KAPPA[k] * gamma[k] + 1.0 / (k as f64 + 1.0);
        }
        Self { gamma, alpha, err_const, ..Self::default() }
    }

    fn change_d(&mut self, factor: f64) {
        let order = self.order;
        let r = compute_r(order, factor);
        let u = compute_r(order, 1.0);
        let ru = r * u;
        let n = self.d[0].len();
        let mut updated = vec![vec![0.0; n]; order + 1];
        for (i, row) in updated.iter_mut().enumerate() {
            for j in 0..=order {
                let c = ru[(j, i)];
                if c != 0.0 {
                    for k in 0..n {
                        row[k] += c * self.d[j][k];
                    }
                }
            }
        }
        for (i, row) in updated.into_iter().enumerate() {
            self.d[i] = row;
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn newton(
        &mut self,
        sys: &dyn OdeSystem,
        t_new: f64,
        predict: &[f64],
        c: f64,
        psi: &[f64],
        lu: &nalgebra::LU<f64, nalgebra::Dyn, nalgebra::Dyn>,
        scale: &[f64],
        tol_newton: f64,
    ) -> (bool, usize, Vec<f64>, Vec<f64>) {
        let n = predict.len();
        let mass = sys.mass();
        let mut d = vec![0.0; n];
        let mut y = predict.to_vec();
        let mut f = vec![0.0; n];
        let mut old_norm: Option<f64> = None;
        let mut converged = false;
        let mut iters = 0;
        for k in 0..NEWTON_MAXITER {
            iters = k + 1;
            sys.rhs(t_new, &y, &mut f);
            self.stats.rhs_evals += 1;
            if !all_finite(&f) {
                break;
            }
            let rhs: Vec<f64> = (0..n).map(|i| c * f[i] - mass[i] * (psi[i] + d[i])).collect();
            let Some(dy) = lu.solve(&DVector::from_vec(rhs)) else { break };
            let norm = (dy.iter().zip(scale).map(|(a, s)| (a / s) * (a / s)).sum::<f64>() / n as f64).sqrt();
            let rate = old_norm.map(|o| norm / o);
            if let Some(r) = rate
                && (r >= 1.0 || r.powi((NEWTON_MAXITER - k) as i32) / (1.0 - r) * norm > tol_newton)
            {
                break;
            }
            for i in 0..n {
                y[i] += dy[i];
                d[i] += dy[i];
            }
            if norm == 0.0 || rate.is_some_and(|r| r / (1.0 - r) * norm < tol_newton) {
                converged = true;
                break;
            }
            old_norm = Some(norm);
        }
        (converged, iters, y, d)
    }
}

impl Integrator for Bdf {
    fn restart(&mut self, sys: &dyn OdeSystem, t: f64, y: &[f64]) {
        let n = y.len();
        self.t = t;
        self.t_prev = t;
        self.y = y.to_vec();
        self.order = 1;
        self.n_equal = 0;
        self.h = 0.0;
        self.d = vec![vec![0.0; n]; MAX_ORDER + 3];
        self.d[0] = y.to_vec();
        self.jac = None;
        self.jac_current = false;
        let _ = sys;
    }

    fn time(&self) -> f64 {
        self.t
    }

    fn state(&self) -> &[f64] {
        &self.y
    }

    fn previous_time(&self) -> f64 {
        self.t_prev
    }

    fn stats(&self) -> Stats {
        self.stats
    }

    fn step(&mut self, sys: &dyn OdeSystem, t_limit: f64, tol: &Tolerances) -> Result<()> {
        let n = self.y.len();
        if n == 0 {
            self.t_prev = self.t;
            self.t = t_limit;
            return Ok(());
        }
        let mass: Vec<f64> = sys.mass().to_vec();
        if self.h <= 0.0 {
            let mut f = vec![0.0; n];
            sys.rhs(self.t, &self.y, &mut f);
            self.stats.rhs_evals += 1;
            self.h = initial_step(sys, self.t, &self.y, &f, tol, 1.0).max(1e-14).min(t_limit - self.t);
            for k in 0..n {
                self.d[1][k] = f[k] * self.h * mass[k];
            }
        }
        if self.h > tol.h_max {
            let factor = tol.h_max / self.h;
            self.change_d(factor);
            self.h = tol.h_max;
            self.n_equal = 0;
        }
        let eps = f64::EPSILON;
        let newton_tol = (10.0 * eps / tol.rtol).max(0.03f64.min(tol.rtol.sqrt()));
        let mut lu: Option<nalgebra::LU<f64, nalgebra::Dyn, nalgebra::Dyn>> = None;
        for _ in 0..tol.max_steps {
            if self.h.abs() < 1e-14 * self.t.abs().max(1e-300) {
                return Err(Error::Integrator(format!("step size underflow at t = {}", self.t)));
            }
            let mut t_new = self.t + self.h;
            // A step that ends within rounding of the limit lands on it, or the next step would be
            // of the order of the rounding error and destroy the difference table. The band stays
            // small, so that a step rejected near the limit can still shrink.
            if t_new >= t_limit - 2e-14 * (t_limit.abs() + self.h.abs()) {
                t_new = t_limit;
                let factor = (t_new - self.t) / self.h;
                self.change_d(factor);
                self.n_equal = 0;
                lu = None;
                self.h = t_new - self.t;
            }
            let h = t_new - self.t;
            let order = self.order;
            let mut predict = vec![0.0; n];
            for j in 0..=order {
                for k in 0..n {
                    predict[k] += self.d[j][k];
                }
            }
            let scale: Vec<f64> = predict.iter().map(|v| tol.atol + tol.rtol * v.abs()).collect();
            let mut psi = vec![0.0; n];
            for j in 1..=order {
                for k in 0..n {
                    psi[k] += self.d[j][k] * self.gamma[j];
                }
            }
            for v in psi.iter_mut() {
                *v /= self.alpha[order];
            }
            let c = h / self.alpha[order];
            let (converged, iters, y_new, d_new) = loop {
                if self.jac.is_none() {
                    let mut j = DMatrix::zeros(n, n);
                    sys.jacobian(t_new, &predict, &mut j);
                    self.stats.jacobians += 1;
                    self.jac = Some(j);
                    self.jac_current = true;
                    lu = None;
                }
                if lu.is_none() {
                    let mut a = -(self.jac.as_ref().expect("jacobian")) * c;
                    for k in 0..n {
                        a[(k, k)] += mass[k];
                    }
                    self.stats.factorizations += 1;
                    lu = Some(a.lu());
                }
                let factor = lu.take().expect("factorization");
                let (ok, its, y, d) = self.newton(sys, t_new, &predict, c, &psi, &factor, &scale, newton_tol);
                lu = Some(factor);
                if ok || self.jac_current {
                    break (ok, its, y, d);
                }
                self.jac = None;
            };
            if !converged {
                self.stats.rejected += 1;
                self.h *= 0.5;
                self.change_d(0.5);
                self.n_equal = 0;
                lu = None;
                continue;
            }
            let safety = 0.9 * (2.0 * NEWTON_MAXITER as f64 + 1.0) / (2.0 * NEWTON_MAXITER as f64 + iters as f64);
            let scale: Vec<f64> = y_new.iter().map(|v| tol.atol + tol.rtol * v.abs()).collect();
            let weighted = |v: &[f64], c: f64| -> f64 {
                (v.iter().zip(&scale).map(|(a, s)| (c * a / s) * (c * a / s)).sum::<f64>() / n as f64).sqrt()
            };
            let err_norm = weighted(&d_new, self.err_const[order]);
            if err_norm > 1.0 {
                self.stats.rejected += 1;
                let factor = (safety * err_norm.powf(-1.0 / (order as f64 + 1.0))).max(0.2);
                self.h *= factor;
                self.change_d(factor);
                self.n_equal = 0;
                lu = None;
                continue;
            }
            // Accepted.
            self.stats.steps += 1;
            self.jac_current = false;
            self.t_prev = self.t;
            self.t = t_new;
            self.y = y_new;
            self.n_equal += 1;
            for k in 0..n {
                self.d[order + 2][k] = d_new[k] - self.d[order + 1][k];
                self.d[order + 1][k] = d_new[k];
            }
            for j in (0..=order).rev() {
                for k in 0..n {
                    let v = self.d[j + 1][k];
                    self.d[j][k] += v;
                }
            }
            if self.n_equal < order + 1 {
                return Ok(());
            }
            let err_m = if order > 1 { weighted(&self.d[order], self.err_const[order - 1]) } else { f64::INFINITY };
            let err_p =
                if order < MAX_ORDER { weighted(&self.d[order + 2], self.err_const[order + 1]) } else { f64::INFINITY };
            let factors = [
                err_m.max(1e-300).powf(-1.0 / order as f64),
                err_norm.max(1e-300).powf(-1.0 / (order as f64 + 1.0)),
                err_p.max(1e-300).powf(-1.0 / (order as f64 + 2.0)),
            ];
            let mut best = 1;
            for k in [0, 2] {
                if factors[k] > factors[best] {
                    best = k;
                }
            }
            self.order = (order as i64 + best as i64 - 1) as usize;
            let factor = (safety * factors[best]).min(10.0);
            self.h *= factor;
            self.change_d(factor);
            self.n_equal = 0;
            return Ok(());
        }
        Err(Error::Integrator(format!("too many step rejections at t = {}", self.t)))
    }

    fn interpolate(&self, t: f64, out: &mut [f64]) {
        // The current backward differences define the interpolating polynomial through the
        // recent solution values, on the current step grid.
        let order = self.order;
        let h = self.h.max(1e-300);
        let n = out.len();
        out.copy_from_slice(&self.d[0]);
        let mut p = 1.0;
        for j in 0..order {
            p *= (t - (self.t - h * j as f64)) / (h * (1.0 + j as f64));
            for k in 0..n {
                out[k] += self.d[j + 1][k] * p;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Stiff linear test problem `y' = -lambda (y - cos t) - sin t` with solution `y = cos t`.
    struct Prothero {
        lambda: f64,
        mass: Vec<f64>,
    }

    impl OdeSystem for Prothero {
        fn dim(&self) -> usize {
            1
        }
        fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
            dy[0] = -self.lambda * (y[0] - t.cos()) - t.sin();
        }
        fn jacobian(&self, _t: f64, _y: &[f64], jac: &mut DMatrix<f64>) {
            jac[(0, 0)] = -self.lambda;
        }
        fn time_derivative(&self, t: f64, _y: &[f64], out: &mut [f64]) {
            out[0] = -self.lambda * t.sin() - t.cos();
        }
        fn mass(&self) -> &[f64] {
            &self.mass
        }
    }

    /// Robertson's chemical kinetics, the classic stiff benchmark, with its conservation law.
    struct Robertson {
        mass: Vec<f64>,
    }

    impl OdeSystem for Robertson {
        fn dim(&self) -> usize {
            3
        }
        fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
            dy[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
            dy[2] = 3e7 * y[1] * y[1];
            dy[1] = -dy[0] - dy[2];
        }
        fn jacobian(&self, _t: f64, y: &[f64], j: &mut DMatrix<f64>) {
            j[(0, 0)] = -0.04;
            j[(0, 1)] = 1e4 * y[2];
            j[(0, 2)] = 1e4 * y[1];
            j[(2, 1)] = 6e7 * y[1];
            j[(1, 0)] = 0.04;
            j[(1, 1)] = -1e4 * y[2] - 6e7 * y[1];
            j[(1, 2)] = -1e4 * y[1];
        }
        fn time_derivative(&self, _t: f64, _y: &[f64], _out: &mut [f64]) {}
        fn mass(&self) -> &[f64] {
            &self.mass
        }
    }

    fn run(int: &mut dyn Integrator, sys: &dyn OdeSystem, y0: &[f64], t_end: f64, tol: &Tolerances) -> Vec<f64> {
        int.restart(sys, 0.0, y0);
        while int.time() < t_end {
            int.step(sys, t_end, tol).unwrap();
        }
        int.state().to_vec()
    }

    #[test]
    fn stiff_linear_problem_reaches_tolerance() {
        let sys = Prothero { lambda: 1e6, mass: vec![1.0] };
        let tol = Tolerances { rtol: 1e-8, atol: 1e-10, ..Tolerances::default() };
        for int in [&mut Rodas4::new() as &mut dyn Integrator, &mut Bdf::new()] {
            let y = run(int, &sys, &[1.0], 5.0, &tol);
            assert!((y[0] - 5f64.cos()).abs() < 1e-6, "{}", y[0]);
        }
    }

    #[test]
    fn robertson_matches_reference_and_conserves_mass() {
        // Reference at t = 40: agrees with the printed cvRoberts_dns output of SUNDIALS
        // (7.158e-1, 9.185e-6, 2.842e-1); docs/VALIDATION.md records the Radau cross-check.
        let reference = [0.7158270688, 9.185534764e-6, 0.2841637457];
        let sys = Robertson { mass: vec![1.0; 3] };
        let tol = Tolerances { rtol: 1e-9, atol: 1e-14, ..Tolerances::default() };
        for int in [&mut Rodas4::new() as &mut dyn Integrator, &mut Bdf::new()] {
            let y = run(int, &sys, &[1.0, 0.0, 0.0], 40.0, &tol);
            for k in 0..3 {
                assert!((y[k] - reference[k]).abs() < 1e-6 * reference[k].max(1e-3), "{k}: {y:?}");
            }
            assert!((y.iter().sum::<f64>() - 1.0).abs() < 1e-10);
        }
    }

    /// Fixed-step order check: the error ratio on halving h gives the method's order.
    #[test]
    fn rodas_converges_with_order_four() {
        struct Nonlinear {
            mass: Vec<f64>,
        }
        impl OdeSystem for Nonlinear {
            fn dim(&self) -> usize {
                2
            }
            fn rhs(&self, t: f64, y: &[f64], dy: &mut [f64]) {
                dy[0] = y[1];
                dy[1] = -y[0] + t.sin() * 0.3;
            }
            fn jacobian(&self, _t: f64, _y: &[f64], j: &mut DMatrix<f64>) {
                j[(0, 1)] = 1.0;
                j[(1, 0)] = -1.0;
            }
            fn time_derivative(&self, t: f64, _y: &[f64], out: &mut [f64]) {
                out[1] = t.cos() * 0.3;
            }
            fn mass(&self) -> &[f64] {
                &self.mass
            }
        }
        let sys = Nonlinear { mass: vec![1.0; 2] };
        let tight = Tolerances { rtol: 1e-13, atol: 1e-14, ..Tolerances::default() };
        let exact = run(&mut Rodas4::new(), &sys, &[1.0, 0.0], 1.0, &tight);
        let mut errors = Vec::new();
        for steps in [8, 16, 32] {
            let h = 1.0 / steps as f64;
            let mut int = Rodas4::new();
            int.restart(&sys, 0.0, &[1.0, 0.0]);
            // Loose tolerances never reject, and the step is capped at the grid.
            let tol = Tolerances { rtol: 1.0, atol: 1.0, h_max: h, ..Tolerances::default() };
            for k in 1..=steps {
                int.h = h;
                int.step(&sys, k as f64 * h, &tol).unwrap();
            }
            let e = int.state().iter().zip(&exact).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
            errors.push(e);
        }
        let order1 = (errors[0] / errors[1]).log2();
        let order2 = (errors[1] / errors[2]).log2();
        assert!(order1 > 3.7 && order2 > 3.7, "orders {order1} {order2} from {errors:?}");
    }

    #[test]
    fn algebraic_equation_is_satisfied() {
        // y0' = -y0, 0 = y1 - 2 y0.
        struct Dae {
            mass: Vec<f64>,
        }
        impl OdeSystem for Dae {
            fn dim(&self) -> usize {
                2
            }
            fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
                dy[0] = -y[0];
                dy[1] = y[1] - 2.0 * y[0];
            }
            fn jacobian(&self, _t: f64, _y: &[f64], j: &mut DMatrix<f64>) {
                j[(0, 0)] = -1.0;
                j[(1, 0)] = -2.0;
                j[(1, 1)] = 1.0;
            }
            fn time_derivative(&self, _t: f64, _y: &[f64], _out: &mut [f64]) {}
            fn mass(&self) -> &[f64] {
                &self.mass
            }
        }
        let sys = Dae { mass: vec![1.0, 0.0] };
        let tol = Tolerances { rtol: 1e-8, atol: 1e-10, ..Tolerances::default() };
        for int in [&mut Rodas4::new() as &mut dyn Integrator, &mut Bdf::new()] {
            let y = run(int, &sys, &[1.0, 2.0], 2.0, &tol);
            assert!((y[0] - (-2f64).exp()).abs() < 1e-6, "{y:?}");
            assert!((y[1] - 2.0 * y[0]).abs() < 1e-7, "{y:?}");
        }
    }
}
