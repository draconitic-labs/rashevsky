//! Solver abstraction for linear, quadratic and mixed-integer programs.
//!
//! A [`Problem`] is a static description. A [`Backend`] opens a [`Session`] on it; a session holds
//! the solver state, so a sequence of bound or objective changes re-solves from the previous basis.
//! Two backends exist: HiGHS (LP, MILP, convex QP; native builds) and [`simplex::Simplex`], a pure
//! Rust dense bounded-variable primal simplex (LP only) used in WebAssembly and as an independent
//! check on HiGHS.

#[cfg(feature = "highs")]
pub mod highs;
pub mod simplex;

use crate::error::Result;
use serde::{Deserialize, Serialize};

/// Direction of optimization.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Sense {
    /// Minimize the objective.
    Minimize,
    /// Maximize the objective.
    Maximize,
}

/// Terminal status of a solve.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Status {
    /// An optimal solution was found.
    Optimal,
    /// The constraints admit no solution.
    Infeasible,
    /// The objective is unbounded over the feasible set.
    Unbounded,
    /// The solver stopped without a proof either way.
    Unknown,
}

/// Mathematical program: `min/max c'x + 0.5 x'Qx` subject to `row_lb <= A x <= row_ub`,
/// `col_lb <= x <= col_ub`, and integrality of the marked columns.
#[derive(Clone, Debug, Default)]
pub struct Problem {
    /// Optimization direction. `None` is treated as minimize.
    pub sense: Option<Sense>,
    /// Linear objective coefficients, one per column.
    pub obj: Vec<f64>,
    /// Column lower bounds.
    pub col_lb: Vec<f64>,
    /// Column upper bounds.
    pub col_ub: Vec<f64>,
    /// Row lower bounds.
    pub row_lb: Vec<f64>,
    /// Row upper bounds.
    pub row_ub: Vec<f64>,
    /// Sparse rows as `(column, coefficient)`.
    pub rows: Vec<Vec<(usize, f64)>>,
    /// Integrality flag per column. Empty means all continuous.
    pub integer: Vec<bool>,
    /// Lower triangle of `Q` as `(row, col, value)` with `row >= col`. Empty means a linear program.
    pub hessian: Vec<(usize, usize, f64)>,
}

impl Problem {
    /// Number of columns.
    pub fn n_cols(&self) -> usize {
        self.obj.len()
    }

    /// Number of rows.
    pub fn n_rows(&self) -> usize {
        self.rows.len()
    }

    /// Adds a column and returns its index.
    pub fn add_col(&mut self, obj: f64, lb: f64, ub: f64, integer: bool) -> usize {
        self.obj.push(obj);
        self.col_lb.push(lb);
        self.col_ub.push(ub);
        if integer || !self.integer.is_empty() {
            self.integer.resize(self.obj.len() - 1, false);
            self.integer.push(integer);
        }
        self.obj.len() - 1
    }

    /// Adds a row and returns its index.
    pub fn add_row(&mut self, lb: f64, ub: f64, coefs: Vec<(usize, f64)>) -> usize {
        self.rows.push(coefs);
        self.row_lb.push(lb);
        self.row_ub.push(ub);
        self.rows.len() - 1
    }

    /// True if any column is integer.
    pub fn is_mip(&self) -> bool {
        self.integer.iter().any(|&b| b)
    }

    /// True if the objective has a quadratic term.
    pub fn is_qp(&self) -> bool {
        !self.hessian.is_empty()
    }

    /// Value of the objective at `x`, including the quadratic term.
    pub fn objective_at(&self, x: &[f64]) -> f64 {
        let linear: f64 = self.obj.iter().zip(x).map(|(c, v)| c * v).sum();
        let mut quad = 0.0;
        for &(i, j, q) in &self.hessian {
            let term = q * x[i] * x[j];
            quad += if i == j { 0.5 * term } else { term };
        }
        linear + quad
    }

    /// Compressed sparse column form `(start, row index, value)`.
    pub fn to_csc(&self) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
        let n = self.n_cols();
        let mut counts = vec![0usize; n + 1];
        for row in &self.rows {
            for &(j, _) in row {
                counts[j + 1] += 1;
            }
        }
        for j in 0..n {
            counts[j + 1] += counts[j];
        }
        let nnz = counts[n];
        let start = counts.clone();
        let mut next = counts;
        let mut index = vec![0usize; nnz];
        let mut value = vec![0f64; nnz];
        for (i, row) in self.rows.iter().enumerate() {
            for &(j, v) in row {
                let k = next[j];
                index[k] = i;
                value[k] = v;
                next[j] += 1;
            }
        }
        (start, index, value)
    }
}

/// Result of a solve.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Solution {
    /// Terminal status.
    pub status: Status,
    /// Objective value in the problem's own sense, or NaN when not optimal.
    pub objective: f64,
    /// Primal column values.
    pub x: Vec<f64>,
    /// Row activities `A x`.
    pub row_activity: Vec<f64>,
    /// Row duals (shadow prices). Empty for mixed-integer solves.
    pub row_dual: Vec<f64>,
    /// Column duals (reduced costs). Empty for mixed-integer solves.
    pub col_dual: Vec<f64>,
}

impl Solution {
    /// True when the status is [`Status::Optimal`].
    pub fn is_optimal(&self) -> bool {
        self.status == Status::Optimal
    }
}

/// Solver state on one problem; changes re-solve from the previous basis where the backend can.
pub trait Session: Send {
    /// Changes the bounds of one column.
    fn set_col_bounds(&mut self, col: usize, lb: f64, ub: f64) -> Result<()>;
    /// Changes the bounds of one row.
    fn set_row_bounds(&mut self, row: usize, lb: f64, ub: f64) -> Result<()>;
    /// Replaces the objective.
    fn set_objective(&mut self, sense: Sense, obj: &[f64]) -> Result<()>;
    /// Appends a row and returns its index.
    fn add_row(&mut self, lb: f64, ub: f64, coefs: &[(usize, f64)]) -> Result<usize>;
    /// Solves the current problem.
    fn solve(&mut self) -> Result<Solution>;
}

/// Factory for solver sessions.
pub trait Backend: Send + Sync {
    /// Name reported in output files.
    fn name(&self) -> &'static str;
    /// Opens a session on `problem`.
    fn open(&self, problem: &Problem) -> Result<Box<dyn Session>>;
}

/// The default backend of this build: HiGHS when compiled in, otherwise the pure Rust simplex.
pub fn default_backend() -> Box<dyn Backend> {
    #[cfg(feature = "highs")]
    {
        Box::new(highs::Highs::default())
    }
    #[cfg(not(feature = "highs"))]
    {
        Box::new(simplex::Simplex::default())
    }
}
