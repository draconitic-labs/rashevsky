//! Pure Rust bounded-variable primal simplex with a dense basis inverse.
//!
//! The method: every row `i` gets a logical variable `s_i` so that `A x - s = 0`, with the row bounds
//! moved to `s`. A composite phase 1 minimizes the sum of bound violations of basic variables, phase
//! 2 optimizes the true objective. The explicit inverse is updated by Gauss-Jordan pivots and rebuilt
//! by LU every [`REFACTOR_EVERY`] iterations. Dantzig pricing falls back to Bland's rule after a run
//! of degenerate pivots, which guarantees termination on the heavily degenerate LPs of metabolic
//! networks. Dense storage limits this backend to problems of a few thousand rows; it exists for
//! WebAssembly builds, where HiGHS is unavailable, and as an independent check on HiGHS.

use super::{Backend, Problem, Sense, Session, Solution, Status};
use crate::error::{Error, Result};
use nalgebra::DMatrix;

const REFACTOR_EVERY: usize = 64;
const STALL_BEFORE_BLAND: usize = 200;

/// Settings of the dense simplex.
#[derive(Clone, Debug)]
pub struct Simplex {
    /// Bound-violation tolerance on primal variables.
    pub primal_tolerance: f64,
    /// Reduced-cost tolerance.
    pub dual_tolerance: f64,
    /// Smallest pivot accepted.
    pub pivot_tolerance: f64,
    /// Iteration limit per solve.
    pub max_iterations: usize,
}

impl Default for Simplex {
    fn default() -> Self {
        Self {
            primal_tolerance: 1e-9,
            dual_tolerance: 1e-9,
            pivot_tolerance: 1e-8,
            max_iterations: 500_000,
        }
    }
}

impl Backend for Simplex {
    fn name(&self) -> &'static str {
        "rashevsky-simplex"
    }

    fn open(&self, problem: &Problem) -> Result<Box<dyn Session>> {
        if problem.is_mip() || problem.is_qp() {
            return Err(Error::Solver(
                "the pure Rust simplex solves linear programs only; build with the `highs` feature".into(),
            ));
        }
        Ok(Box::new(SimplexSession::new(problem, self.clone())))
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum VarStat {
    Basic,
    AtLower,
    AtUpper,
    Free,
}

/// Simplex state. Variables `0..n` are structural, `n..n+m` logical.
pub struct SimplexSession {
    opts: Simplex,
    n: usize,
    m: usize,
    cols: Vec<Vec<(usize, f64)>>,
    lb: Vec<f64>,
    ub: Vec<f64>,
    cost: Vec<f64>,
    negate: bool,
    stat: Vec<VarStat>,
    basis: Vec<usize>,
    binv: DMatrix<f64>,
    x: Vec<f64>,
    dirty_basis: bool,
}

impl SimplexSession {
    fn new(p: &Problem, opts: Simplex) -> Self {
        let n = p.n_cols();
        let m = p.n_rows();
        let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n + m];
        for (i, row) in p.rows.iter().enumerate() {
            for &(j, v) in row {
                cols[j].push((i, v));
            }
        }
        for i in 0..m {
            cols[n + i].push((i, -1.0));
        }
        let mut lb = p.col_lb.clone();
        let mut ub = p.col_ub.clone();
        lb.extend_from_slice(&p.row_lb);
        ub.extend_from_slice(&p.row_ub);
        let negate = p.sense == Some(Sense::Maximize);
        let mut cost: Vec<f64> = p.obj.iter().map(|&c| if negate { -c } else { c }).collect();
        cost.resize(n + m, 0.0);
        let mut session = Self {
            opts,
            n,
            m,
            cols,
            lb,
            ub,
            cost,
            negate,
            stat: vec![VarStat::AtLower; n + m],
            basis: (n..n + m).collect(),
            binv: DMatrix::identity(m, m),
            x: vec![0.0; n + m],
            dirty_basis: false,
        };
        session.reset_basis();
        session
    }

    fn reset_basis(&mut self) {
        for j in 0..self.n {
            self.place_nonbasic(j);
        }
        for i in 0..self.m {
            self.stat[self.n + i] = VarStat::Basic;
            self.basis[i] = self.n + i;
        }
        self.binv = -DMatrix::<f64>::identity(self.m, self.m);
        self.dirty_basis = false;
    }

    /// Puts a nonbasic variable on the finite bound nearest zero, or at zero when free.
    fn place_nonbasic(&mut self, j: usize) {
        let (l, u) = (self.lb[j], self.ub[j]);
        let (stat, value) = match (l.is_finite(), u.is_finite()) {
            (true, true) => {
                if l.abs() <= u.abs() {
                    (VarStat::AtLower, l)
                } else {
                    (VarStat::AtUpper, u)
                }
            }
            (true, false) => (VarStat::AtLower, l),
            (false, true) => (VarStat::AtUpper, u),
            (false, false) => (VarStat::Free, 0.0),
        };
        self.stat[j] = stat;
        self.x[j] = value;
    }

    /// Rebuilds the inverse from scratch. Returns false if the basis is singular.
    fn refactor(&mut self) -> bool {
        let mut b = DMatrix::<f64>::zeros(self.m, self.m);
        for (k, &var) in self.basis.iter().enumerate() {
            for &(i, v) in &self.cols[var] {
                b[(i, k)] = v;
            }
        }
        match b.lu().try_inverse() {
            Some(inv) => {
                self.binv = inv;
                true
            }
            None => false,
        }
    }

    /// Recomputes the basic values from the nonbasic ones.
    fn compute_basics(&mut self) {
        let mut rhs = vec![0.0; self.m];
        for j in 0..self.n + self.m {
            if self.stat[j] != VarStat::Basic && self.x[j] != 0.0 {
                for &(i, v) in &self.cols[j] {
                    rhs[i] -= v * self.x[j];
                }
            }
        }
        for k in 0..self.m {
            let mut value = 0.0;
            for i in 0..self.m {
                value += self.binv[(k, i)] * rhs[i];
            }
            self.x[self.basis[k]] = value;
        }
    }

    fn column_ftran(&self, j: usize) -> Vec<f64> {
        let mut alpha = vec![0.0; self.m];
        for &(i, v) in &self.cols[j] {
            for (k, a) in alpha.iter_mut().enumerate() {
                *a += self.binv[(k, i)] * v;
            }
        }
        alpha
    }

    fn infeasibility(&self) -> f64 {
        let tol = self.opts.primal_tolerance;
        let mut total = 0.0;
        for &var in &self.basis {
            let v = self.x[var];
            if v < self.lb[var] - tol {
                total += self.lb[var] - v;
            } else if v > self.ub[var] + tol {
                total += v - self.ub[var];
            }
        }
        total
    }

    /// Runs the simplex. The status is that of the final phase.
    fn run(&mut self) -> Status {
        if self.dirty_basis {
            for j in 0..self.n + self.m {
                if self.stat[j] != VarStat::Basic {
                    self.settle_nonbasic(j);
                }
            }
            if !self.refactor() {
                self.reset_basis();
            }
            self.dirty_basis = false;
        }
        self.compute_basics();
        let tol = self.opts.primal_tolerance;
        let dtol = self.opts.dual_tolerance;
        let mut since_refactor = 0;
        let mut degenerate = 0usize;
        for _ in 0..self.opts.max_iterations {
            if since_refactor >= REFACTOR_EVERY {
                if !self.refactor() {
                    self.reset_basis();
                }
                self.compute_basics();
                since_refactor = 0;
            }
            let phase1 = self.infeasibility() > 0.0;
            // Basic cost vector and simplex multipliers.
            let mut y = vec![0.0; self.m];
            for k in 0..self.m {
                let var = self.basis[k];
                let c = if phase1 {
                    let v = self.x[var];
                    if v < self.lb[var] - tol {
                        -1.0
                    } else if v > self.ub[var] + tol {
                        1.0
                    } else {
                        0.0
                    }
                } else {
                    self.cost[var]
                };
                if c != 0.0 {
                    for (i, yi) in y.iter_mut().enumerate() {
                        *yi += c * self.binv[(k, i)];
                    }
                }
            }
            // Pricing.
            let bland = degenerate >= STALL_BEFORE_BLAND;
            let mut entering: Option<(usize, f64)> = None;
            let mut best = 0.0;
            for j in 0..self.n + self.m {
                let stat = self.stat[j];
                if stat == VarStat::Basic || self.lb[j] == self.ub[j] {
                    continue;
                }
                let mut d = if phase1 { 0.0 } else { self.cost[j] };
                for &(i, v) in &self.cols[j] {
                    d -= y[i] * v;
                }
                let dir = match stat {
                    VarStat::AtLower if d < -dtol => 1.0,
                    VarStat::AtUpper if d > dtol => -1.0,
                    VarStat::Free if d.abs() > dtol => -d.signum(),
                    _ => continue,
                };
                if bland {
                    entering = Some((j, dir));
                    break;
                }
                if d.abs() > best {
                    best = d.abs();
                    entering = Some((j, dir));
                }
            }
            let Some((enter, dir)) = entering else {
                return if phase1 { Status::Infeasible } else { Status::Optimal };
            };
            let alpha = self.column_ftran(enter);
            // Two-pass ratio test: the largest pivot among the blocking steps within tolerance.
            let mut t_limit = self.ub[enter] - self.lb[enter];
            if !t_limit.is_finite() {
                t_limit = f64::INFINITY;
            }
            let ptol = self.opts.pivot_tolerance;
            let mut t_relaxed = t_limit;
            for k in 0..self.m {
                let rate = -dir * alpha[k];
                if rate.abs() <= ptol {
                    continue;
                }
                let var = self.basis[k];
                if let Some(t) = self.blocking_step(var, rate, tol, phase1) {
                    t_relaxed = t_relaxed.min(t);
                }
            }
            if !t_relaxed.is_finite() {
                return if phase1 { Status::Unknown } else { Status::Unbounded };
            }
            let mut leave: Option<usize> = None;
            let mut pivot_best = 0.0;
            let mut t_step = t_limit;
            for k in 0..self.m {
                let rate = -dir * alpha[k];
                if rate.abs() <= ptol {
                    continue;
                }
                let var = self.basis[k];
                if let Some(t_exact) = self.blocking_step(var, rate, 0.0, phase1)
                    && t_exact <= t_relaxed
                {
                    let better = if bland {
                        match leave {
                            None => true,
                            Some(l) => var < self.basis[l],
                        }
                    } else {
                        alpha[k].abs() > pivot_best
                    };
                    if better {
                        pivot_best = alpha[k].abs();
                        leave = Some(k);
                        t_step = t_exact.max(0.0);
                    }
                }
            }
            if leave.is_some() && t_limit <= t_step {
                leave = None;
                t_step = t_limit;
            }
            if t_step <= 1e-12 {
                degenerate += 1;
            } else {
                degenerate = 0;
            }
            // Apply the step.
            for k in 0..self.m {
                let var = self.basis[k];
                self.x[var] += -dir * alpha[k] * t_step;
            }
            self.x[enter] += dir * t_step;
            match leave {
                None => {
                    // The entering variable moved to its opposite bound.
                    let (stat, value) = if dir > 0.0 {
                        (VarStat::AtUpper, self.ub[enter])
                    } else {
                        (VarStat::AtLower, self.lb[enter])
                    };
                    self.stat[enter] = stat;
                    self.x[enter] = value;
                }
                Some(r) => {
                    let out = self.basis[r];
                    self.settle_leaving(out, -dir * alpha[r]);
                    self.basis[r] = enter;
                    self.stat[enter] = VarStat::Basic;
                    self.pivot(r, &alpha);
                    since_refactor += 1;
                }
            }
        }
        Status::Unknown
    }

    /// Step length at which basic `var`, moving at `rate`, reaches its blocking bound.
    fn blocking_step(&self, var: usize, rate: f64, tol: f64, phase1: bool) -> Option<f64> {
        let v = self.x[var];
        let (l, u) = (self.lb[var], self.ub[var]);
        let ftol = self.opts.primal_tolerance;
        if phase1 {
            // An infeasible variable blocks at the violated bound once it has moved back to it.
            if v < l - ftol {
                return if rate > 0.0 { Some((l - v + tol) / rate) } else { None };
            }
            if v > u + ftol {
                return if rate < 0.0 { Some((u - v - tol) / rate) } else { None };
            }
        }
        if rate > 0.0 {
            if u.is_finite() { Some(((u - v + tol) / rate).max(0.0)) } else { None }
        } else if l.is_finite() {
            Some(((l - v - tol) / rate).max(0.0))
        } else {
            None
        }
    }

    fn settle_leaving(&mut self, var: usize, rate: f64) {
        let (l, u) = (self.lb[var], self.ub[var]);
        let v = self.x[var];
        if rate > 0.0 {
            if u.is_finite() && (v >= u - 1e-7 || !l.is_finite()) {
                self.stat[var] = VarStat::AtUpper;
                self.x[var] = u;
            } else if l.is_finite() {
                self.stat[var] = VarStat::AtLower;
                self.x[var] = l;
            } else {
                self.stat[var] = VarStat::Free;
            }
        } else if l.is_finite() && (v <= l + 1e-7 || !u.is_finite()) {
            self.stat[var] = VarStat::AtLower;
            self.x[var] = l;
        } else if u.is_finite() {
            self.stat[var] = VarStat::AtUpper;
            self.x[var] = u;
        } else {
            self.stat[var] = VarStat::Free;
        }
    }

    fn pivot(&mut self, r: usize, alpha: &[f64]) {
        let m = self.m;
        let p = alpha[r];
        for i in 0..m {
            self.binv[(r, i)] /= p;
        }
        for k in 0..m {
            if k != r && alpha[k] != 0.0 {
                let factor = alpha[k];
                for i in 0..m {
                    let delta = factor * self.binv[(r, i)];
                    self.binv[(k, i)] -= delta;
                }
            }
        }
    }

    /// Moves a nonbasic variable onto a bound that exists after a bound change.
    fn settle_nonbasic(&mut self, j: usize) {
        let (l, u) = (self.lb[j], self.ub[j]);
        match self.stat[j] {
            VarStat::AtLower if l.is_finite() => self.x[j] = l,
            VarStat::AtUpper if u.is_finite() => self.x[j] = u,
            _ => self.place_nonbasic(j),
        }
    }

    fn extract(&self, status: Status) -> Solution {
        let mut row_activity = vec![0.0; self.m];
        for (j, value) in self.x.iter().take(self.n).enumerate() {
            for &(i, v) in &self.cols[j] {
                row_activity[i] += v * value;
            }
        }
        let sign = if self.negate { -1.0 } else { 1.0 };
        let objective: f64 = self
            .x
            .iter()
            .take(self.n)
            .zip(&self.cost)
            .map(|(x, c)| x * c)
            .sum::<f64>()
            * sign;
        let mut y = vec![0.0; self.m];
        for k in 0..self.m {
            let c = self.cost[self.basis[k]];
            if c != 0.0 {
                for (i, yi) in y.iter_mut().enumerate() {
                    *yi += c * self.binv[(k, i)];
                }
            }
        }
        let col_dual: Vec<f64> = (0..self.n)
            .map(|j| {
                let mut d = self.cost[j];
                for &(i, v) in &self.cols[j] {
                    d -= y[i] * v;
                }
                d * sign
            })
            .collect();
        let row_dual: Vec<f64> = y.iter().map(|v| v * sign).collect();
        Solution {
            status,
            objective: if status == Status::Optimal { objective } else { f64::NAN },
            x: self.x[..self.n].to_vec(),
            row_activity,
            row_dual,
            col_dual,
        }
    }
}

impl Session for SimplexSession {
    fn set_col_bounds(&mut self, col: usize, lb: f64, ub: f64) -> Result<()> {
        self.lb[col] = lb;
        self.ub[col] = ub;
        self.dirty_basis = true;
        Ok(())
    }

    fn set_row_bounds(&mut self, row: usize, lb: f64, ub: f64) -> Result<()> {
        self.lb[self.n + row] = lb;
        self.ub[self.n + row] = ub;
        self.dirty_basis = true;
        Ok(())
    }

    fn set_objective(&mut self, sense: Sense, obj: &[f64]) -> Result<()> {
        if obj.len() != self.n {
            return Err(Error::Solver("objective length differs from column count".into()));
        }
        self.negate = sense == Sense::Maximize;
        for (c, &v) in self.cost.iter_mut().zip(obj) {
            *c = if self.negate { -v } else { v };
        }
        self.dirty_basis = true;
        Ok(())
    }

    fn add_row(&mut self, lb: f64, ub: f64, coefs: &[(usize, f64)]) -> Result<usize> {
        // Logical variables sit behind the structural ones, so insert at the end of the logicals.
        let row = self.m;
        self.m += 1;
        for &(j, v) in coefs {
            self.cols[j].push((row, v));
        }
        self.cols.push(vec![(row, -1.0)]);
        self.lb.push(lb);
        self.ub.push(ub);
        self.cost.push(0.0);
        self.stat.push(VarStat::Basic);
        self.x.push(0.0);
        self.basis.push(self.n + row);
        let old = std::mem::replace(&mut self.binv, DMatrix::zeros(0, 0));
        self.binv = old.insert_row(row, 0.0).insert_column(row, 0.0);
        self.dirty_basis = true;
        Ok(row)
    }

    fn solve(&mut self) -> Result<Solution> {
        for j in 0..self.n + self.m {
            if self.lb[j] > self.ub[j] {
                return Ok(self.extract(Status::Infeasible));
            }
        }
        let status = self.run();
        Ok(self.extract(status))
    }
}
