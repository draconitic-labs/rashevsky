//! Steady states of kinetic models with moiety conservation.
//!
//! The stoichiometric matrix `N` factors as `N = L Nr`, with `Nr` the independent rows and `L` the
//! link matrix (Reder 1988). Dependent species follow from the independent ones through the
//! conserved totals, and damped Newton solves `Nr v(y(s)) = 0` over the independent amounts `s`.

use super::compile::{Coef, KineticModel};
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
    /// Integrate over this time first, which selects the attractor reached from the initial state.
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
}

/// Solves for a steady state.
pub fn steady_state(model: &KineticModel, options: &SteadyOptions) -> Result<SteadyState> {
    let mut y = model.initial_state();
    if options.integrate_first > 0.0 {
        let run = simulate(
            model,
            &[0.0, options.integrate_first],
            &SimOptions { method: Method::Rodas4, rtol: 1e-8, atol: 1e-10, ..SimOptions::default() },
        )?;
        let vals = &run.values[1];
        y = model.state_slots.iter().map(|&s| vals[s]).collect();
    }
    let sys = KineticSystem::new(model, false);
    let m = moieties(model)?;
    let n = y.len();
    let y0 = y.clone();
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
        sys.rhs(0.0, &full, &mut f);
        let norm = f.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
        (f, norm)
    };
    let mut s: Vec<f64> = ind.iter().map(|&i| y[i]).collect();
    let (mut f, mut norm) = residual_of(&s);
    let mut iterations = 0;
    // Newton runs until it stops improving, not merely to the tolerance: control coefficients divide
    // by small fluxes, so they need the steady state as exact as the arithmetic allows.
    while norm > 0.0 && iterations < options.max_iterations {
        iterations += 1;
        let full = lift(&s);
        let mut jac = DMatrix::zeros(n, n);
        sys.jacobian(0.0, &full, &mut jac);
        let mut reduced = DMatrix::zeros(ind.len(), ind.len());
        let jl = &jac * &m.link;
        for (row, &i) in ind.iter().enumerate() {
            reduced.row_mut(row).copy_from(&jl.row(i));
        }
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
    let values = sys.values(0.0, &full);
    let fluxes = model.reactions.iter().map(|r| values[r.rate]).collect();
    Ok(SteadyState { y: full, values, fluxes, residual: norm, iterations })
}
