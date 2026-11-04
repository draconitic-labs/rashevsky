//! Uniform sampling of the steady-state flux polytope, with convergence diagnostics.
//!
//! [`chrr`] is coordinate hit-and-run with rounding (Haraldsdottir et al. 2017): reactions that cannot
//! vary are fixed, the polytope is written `A z <= b` over a null-space basis, and hit-and-run moves
//! along coordinates of a rounded basis. Rounding only affects mixing, since the target stays uniform.
//! This implementation rounds by iterating to isotropic position from sampled covariances (Lovasz and
//! Vempala 2006) rather than with the maximum-volume ellipsoid used by the original authors.
//!
//! [`optgp`] follows Megchelenbrink, Huynen and Marchiori (2014): chains of artificial-centering
//! hit-and-run that all start at the centre of FVA-derived warmup points, with periodic projection
//! back onto the null space. As artificial centering adapts its directions, it is a heuristic sampler
//! whose output should be checked with the diagnostics here.

use super::fva::{FvaOptions, fva_with};
use crate::error::{Error, Result};
use crate::linalg::rref;
use crate::lp::{self, Problem, Sense};
use crate::model::Model;
use nalgebra::{DMatrix, DVector};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};

/// Settings of the samplers.
#[derive(Clone, Copy, Debug)]
pub struct SamplerOptions {
    /// Number of independent chains.
    pub chains: usize,
    /// Stored samples per chain.
    pub samples_per_chain: usize,
    /// Moves between stored samples.
    pub thinning: usize,
    /// Moves discarded at the start of each chain.
    pub burn_in: usize,
    /// Seed of the random streams.
    pub seed: u64,
    /// Fraction of the optimal objective the samples must retain; `None` leaves the objective free.
    pub fraction_of_optimum: Option<f64>,
}

impl Default for SamplerOptions {
    fn default() -> Self {
        Self { chains: 4, samples_per_chain: 1000, thinning: 100, burn_in: 2000, seed: 1, fraction_of_optimum: None }
    }
}

/// Samples and their diagnostics.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FluxSamples {
    /// Reaction identifiers.
    pub reactions: Vec<String>,
    /// `chains[c][s][r]`: reaction `r` in sample `s` of chain `c`.
    pub chains: Vec<Vec<Vec<f64>>>,
    /// Convergence diagnostics.
    pub diagnostics: Diagnostics,
    /// Dimension of the polytope after removing fixed reactions.
    pub dimension: usize,
}

/// Convergence diagnostics over reactions that vary.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Diagnostics {
    /// Largest split R-hat (Gelman et al. 2013).
    pub max_split_rhat: f64,
    /// Smallest effective sample size over all chains.
    pub min_ess: f64,
    /// Largest absolute Geweke z-score, first 10 % against last 50 % of each chain.
    pub max_geweke_z: f64,
    /// Sample mean per reaction.
    pub mean: Vec<f64>,
    /// Sample standard deviation per reaction.
    pub sd: Vec<f64>,
    /// Monte Carlo standard error of the mean per reaction.
    pub mcse: Vec<f64>,
}

fn mean_var(x: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let m = x.iter().sum::<f64>() / n;
    let v = x.iter().map(|a| (a - m) * (a - m)).sum::<f64>() / (n - 1.0).max(1.0);
    (m, v)
}

/// Split R-hat of one quantity across chains.
pub fn split_rhat(chains: &[Vec<f64>]) -> f64 {
    let mut parts: Vec<&[f64]> = Vec::with_capacity(chains.len() * 2);
    for c in chains {
        let half = c.len() / 2;
        parts.push(&c[..half]);
        parts.push(&c[half..2 * half]);
    }
    let n = parts[0].len() as f64;
    if n < 2.0 {
        return f64::NAN;
    }
    let stats: Vec<(f64, f64)> = parts.iter().map(|p| mean_var(p)).collect();
    let within = stats.iter().map(|s| s.1).sum::<f64>() / stats.len() as f64;
    let (_, between_means) = mean_var(&stats.iter().map(|s| s.0).collect::<Vec<_>>());
    if within <= 1e-300 {
        return 1.0;
    }
    (((n - 1.0) / n * within + between_means) / within).sqrt()
}

/// Effective sample size from autocorrelations, summed over Geyer's initial positive sequence.
pub fn effective_sample_size(chains: &[Vec<f64>]) -> f64 {
    let n = chains[0].len();
    let m = chains.len() as f64;
    let stats: Vec<(f64, f64)> = chains.iter().map(|c| mean_var(c)).collect();
    let within = stats.iter().map(|s| s.1).sum::<f64>() / m;
    let (_, between_means) = mean_var(&stats.iter().map(|s| s.0).collect::<Vec<_>>());
    let var_plus = (n as f64 - 1.0) / n as f64 * within + if chains.len() > 1 { between_means } else { 0.0 };
    if var_plus <= 1e-300 {
        return m * n as f64;
    }
    let rho = |lag: usize| -> f64 {
        let mut gamma = 0.0;
        for (c, (mean, _)) in chains.iter().zip(&stats) {
            let mut acc = 0.0;
            for k in 0..n - lag {
                acc += (c[k] - mean) * (c[k + lag] - mean);
            }
            gamma += acc / n as f64;
        }
        1.0 - (within - gamma / m) / var_plus
    };
    let mut sum = 0.0;
    let mut lag = 1;
    while lag + 1 < n {
        let pair = rho(lag) + rho(lag + 1);
        if pair < 0.0 {
            break;
        }
        sum += pair;
        lag += 2;
    }
    (m * n as f64 / (1.0 + 2.0 * sum)).min(m * n as f64)
}

/// Geweke z-score of one chain.
pub fn geweke(chain: &[f64]) -> f64 {
    let n = chain.len();
    let (a, _) = (&chain[..n / 10], ());
    let b = &chain[n / 2..];
    if a.len() < 2 || b.len() < 2 {
        return 0.0;
    }
    let spectral = |x: &[f64]| -> f64 {
        // Variance of the mean corrected by a Bartlett window on the autocovariances.
        let (m, v) = mean_var(x);
        let k = x.len();
        let lags = ((k as f64).sqrt() as usize).max(1);
        let mut s = v;
        for lag in 1..=lags {
            let mut acc = 0.0;
            for i in 0..k - lag {
                acc += (x[i] - m) * (x[i + lag] - m);
            }
            s += 2.0 * (1.0 - lag as f64 / (lags as f64 + 1.0)) * acc / k as f64;
        }
        s.max(0.0) / k as f64
    };
    let (ma, _) = mean_var(a);
    let (mb, _) = mean_var(b);
    let denom = (spectral(a) + spectral(b)).sqrt();
    if denom <= 1e-300 { 0.0 } else { (ma - mb) / denom }
}

fn diagnose(chains: &[Vec<Vec<f64>>], n_reactions: usize) -> Diagnostics {
    let mut max_rhat = 1.0f64;
    let mut min_ess = f64::INFINITY;
    let mut max_z = 0.0f64;
    let mut mean = vec![0.0; n_reactions];
    let mut sd = vec![0.0; n_reactions];
    let mut mcse = vec![0.0; n_reactions];
    for r in 0..n_reactions {
        let series: Vec<Vec<f64>> = chains.iter().map(|c| c.iter().map(|s| s[r]).collect()).collect();
        let all: Vec<f64> = series.iter().flatten().copied().collect();
        let (m, v) = mean_var(&all);
        mean[r] = m;
        sd[r] = v.sqrt();
        if v <= 1e-16 * (1.0 + m * m) {
            continue;
        }
        let ess = effective_sample_size(&series);
        mcse[r] = (v / ess.max(1.0)).sqrt();
        min_ess = min_ess.min(ess);
        let rhat = split_rhat(&series);
        if rhat.is_finite() {
            max_rhat = max_rhat.max(rhat);
        }
        for s in &series {
            max_z = max_z.max(geweke(s).abs());
        }
    }
    Diagnostics {
        max_split_rhat: max_rhat,
        min_ess: if min_ess.is_finite() { min_ess } else { 0.0 },
        max_geweke_z: max_z,
        mean,
        sd,
        mcse,
    }
}

/// Bounds and stoichiometry with the objective restriction applied, as bounds on an extra row.
struct Polytope {
    /// Reactions free to vary.
    free: Vec<usize>,
    /// Value of each reaction fixed by its range.
    fixed: Vec<f64>,
    lb: Vec<f64>,
    ub: Vec<f64>,
    /// `equalities x free`, right-hand side from the fixed reactions.
    eq: DMatrix<f64>,
    eq_rhs: DVector<f64>,
    /// Extra inequality on the objective: `coefs . v >= value`, when requested.
    objective_row: Option<(Vec<f64>, f64)>,
}

fn build_polytope(model: &Model, fraction: Option<f64>, backend: &dyn lp::Backend) -> Result<Polytope> {
    let ranges = fva_with(
        model,
        &FvaOptions { fraction_of_optimum: fraction.unwrap_or(0.0), reactions: None },
        backend,
    )?;
    let n = model.reactions.len();
    let mut fixed = vec![f64::NAN; n];
    let mut free = Vec::new();
    let mut lb = Vec::new();
    let mut ub = Vec::new();
    for r in &ranges {
        if !r.minimum.is_finite() || !r.maximum.is_finite() {
            return Err(Error::Invalid(format!("reaction {} is unbounded; sampling needs a bounded polytope", r.id)));
        }
        if r.maximum - r.minimum < 1e-9 {
            fixed[r.reaction] = 0.5 * (r.minimum + r.maximum);
        } else {
            free.push(r.reaction);
            // Bounds from FVA are tight, which helps rounding, and they imply the model's own.
            lb.push(r.minimum);
            ub.push(r.maximum);
        }
    }
    let m = model.metabolites.len();
    let mut eq = DMatrix::zeros(m, free.len());
    let mut rhs = DVector::zeros(m);
    let mut column_of = vec![usize::MAX; n];
    for (k, &j) in free.iter().enumerate() {
        column_of[j] = k;
    }
    for (j, r) in model.reactions.iter().enumerate() {
        for &(i, c) in &r.stoichiometry {
            if model.metabolites[i].boundary {
                continue;
            }
            if column_of[j] != usize::MAX {
                eq[(i, column_of[j])] = c;
            } else {
                rhs[i] -= c * fixed[j];
            }
        }
    }
    // Keep only independent equations.
    let reduced = rref(&eq.transpose(), 1e-10);
    let mut rows = reduced.pivot_cols.clone();
    rows.sort_unstable();
    let eq_ind = DMatrix::from_fn(rows.len(), free.len(), |a, b| eq[(rows[a], b)]);
    let rhs_ind = DVector::from_iterator(rows.len(), rows.iter().map(|&a| rhs[a]));
    let objective_row = fraction.map(|f| {
        let lp_problem = model.to_lp();
        let mut c = vec![0.0; free.len()];
        for (k, &j) in free.iter().enumerate() {
            c[k] = lp_problem.obj[j];
        }
        (c, f)
    });
    Ok(Polytope { free, fixed, lb, ub, eq: eq_ind, eq_rhs: rhs_ind, objective_row })
}

/// Orthonormal basis of the null space by singular value decomposition.
fn orthonormal_nullspace(a: &DMatrix<f64>) -> DMatrix<f64> {
    let n = a.ncols();
    if a.nrows() == 0 {
        return DMatrix::identity(n, n);
    }
    // The null space of `A` is the eigenspace of zero for `A' A`.
    let ata = a.transpose() * a;
    let eig = ata.symmetric_eigen();
    let scale = eig.eigenvalues.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-300);
    let keep: Vec<usize> = (0..n).filter(|&k| eig.eigenvalues[k].abs() < 1e-11 * scale).collect();
    DMatrix::from_fn(n, keep.len(), |i, k| eig.eigenvectors[(i, keep[k])])
}

/// A point in the interior: the centre of the largest inscribed ball, from a linear program.
fn chebyshev_center(pt: &Polytope, null: &DMatrix<f64>, x0: &DVector<f64>, backend: &dyn lp::Backend) -> Result<DVector<f64>> {
    let d = null.ncols();
    let mut p = Problem { sense: Some(Sense::Maximize), ..Problem::default() };
    for _ in 0..d {
        p.add_col(0.0, f64::NEG_INFINITY, f64::INFINITY, false);
    }
    let radius = p.add_col(1.0, 0.0, f64::INFINITY, false);
    for k in 0..pt.free.len() {
        let row = null.row(k);
        let norm = row.norm();
        let mut coefs: Vec<(usize, f64)> = (0..d).filter(|&c| row[c].abs() > 1e-13).map(|c| (c, row[c])).collect();
        coefs.push((radius, norm));
        p.add_row(f64::NEG_INFINITY, pt.ub[k] - x0[k], coefs.clone());
        let mut lower: Vec<(usize, f64)> = (0..d).filter(|&c| row[c].abs() > 1e-13).map(|c| (c, -row[c])).collect();
        lower.push((radius, norm));
        p.add_row(f64::NEG_INFINITY, x0[k] - pt.lb[k], lower);
    }
    let mut session = backend.open(&p)?;
    let s = session.solve()?;
    if !s.is_optimal() || s.x[radius] <= 0.0 {
        return Err(Error::Invalid("the flux polytope has no interior after removing fixed reactions".into()));
    }
    Ok(DVector::from_column_slice(&s.x[..d]))
}

fn lift(pt: &Polytope, n: usize, v_free: &DVector<f64>) -> Vec<f64> {
    let mut out = pt.fixed.clone();
    debug_assert_eq!(out.len(), n);
    for (k, &j) in pt.free.iter().enumerate() {
        out[j] = v_free[k];
    }
    out
}

/// Coordinate hit-and-run with rounding.
pub fn chrr(model: &Model, options: &SamplerOptions) -> Result<FluxSamples> {
    let backend = lp::default_backend();
    let pt = build_polytope(model, options.fraction_of_optimum, backend.as_ref())?;
    if pt.free.is_empty() {
        return Err(Error::Invalid("no reaction can vary; there is nothing to sample".into()));
    }
    // A particular solution of the equalities by least squares.
    let eq = &pt.eq;
    let x0: DVector<f64> = if eq.nrows() == 0 {
        DVector::zeros(pt.free.len())
    } else {
        let gram = eq * eq.transpose();
        let y = gram.lu().solve(&pt.eq_rhs).ok_or_else(|| Error::Solver("singular equality system".into()))?;
        eq.transpose() * y
    };
    let mut basis = orthonormal_nullspace(eq);
    let d = basis.ncols();
    if d == 0 {
        return Err(Error::Invalid("the flux polytope has dimension zero".into()));
    }
    let mut z = chebyshev_center(&pt, &basis, &x0, backend.as_ref())?;
    let objective_coefs = pt.objective_row.as_ref().map(|(c, _)| DVector::from_column_slice(c));
    let objective_floor = match (&pt.objective_row, options.fraction_of_optimum) {
        (Some(_), Some(f)) => {
            let fba = super::fba::fba_with(model, backend.as_ref())?;
            let c = model.objective_vector();
            let fixed_part: f64 = (0..model.reactions.len())
                .filter(|&j| pt.fixed[j].is_finite())
                .map(|j| c[j] * pt.fixed[j])
                .sum();
            Some(f * fba.objective - fixed_part)
        }
        _ => None,
    };

    // Coordinate hit-and-run in the current basis; `v = x0 + basis z`.
    let run = |basis: &DMatrix<f64>,
               start: &DVector<f64>,
               rng: &mut ChaCha8Rng,
               moves: usize,
               thin: usize,
               store: &mut dyn FnMut(&DVector<f64>)| {
        let k_free = pt.free.len();
        let mut v = &x0 + basis * start;
        for step in 1..=moves {
            let c = rng.random_range(0..d);
            let dir = basis.column(c);
            let mut lo = f64::NEG_INFINITY;
            let mut hi = f64::INFINITY;
            for k in 0..k_free {
                let a = dir[k];
                if a.abs() > 1e-13 {
                    let (t1, t2) = ((pt.lb[k] - v[k]) / a, (pt.ub[k] - v[k]) / a);
                    let (tl, th) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
                    lo = lo.max(tl);
                    hi = hi.min(th);
                }
            }
            if let (Some(obj), Some(floor)) = (&objective_coefs, objective_floor) {
                let slope = obj.dot(&dir);
                let value = obj.dot(&v);
                if slope.abs() > 1e-13 {
                    let t = (floor - value) / slope;
                    if slope > 0.0 {
                        lo = lo.max(t);
                    } else {
                        hi = hi.min(t);
                    }
                }
            }
            if hi > lo {
                let t = lo + rng.random::<f64>() * (hi - lo);
                v.axpy(t, &dir, 1.0);
            }
            if step % thin == 0 {
                store(&v);
            }
        }
    };

    // Rounding: move to isotropic position from sampled covariance, repeated until near isotropic.
    let mut rng = ChaCha8Rng::seed_from_u64(options.seed);
    let rounds = 6;
    let round_samples = (8 * d).clamp(200, 4000);
    for _ in 0..rounds {
        let mut cloud: Vec<DVector<f64>> = Vec::with_capacity(round_samples);
        let thin = (4 * d).max(10);
        let mut collect = |v: &DVector<f64>| cloud.push(basis.transpose() * (v - &x0));
        run(&basis, &z, &mut rng, round_samples * thin, thin, &mut collect);
        let n_c = cloud.len() as f64;
        let mean = cloud.iter().fold(DVector::zeros(d), |acc, p| acc + p) / n_c;
        let mut cov = DMatrix::zeros(d, d);
        for p in &cloud {
            let e = p - &mean;
            cov += &e * e.transpose();
        }
        cov /= n_c - 1.0;
        let eig = cov.clone().symmetric_eigen();
        let largest = eig.eigenvalues.iter().fold(0.0f64, |m, v| m.max(*v)).max(1e-300);
        // Directions with variance far below the largest stay scaled to keep the map invertible.
        let scales = eig.eigenvalues.map(|v| v.max(1e-12 * largest).sqrt());
        let transform = &eig.eigenvectors * DMatrix::from_diagonal(&scales);
        let new_basis = &basis * &transform;
        let inverse = transform.clone().lu().try_inverse().ok_or_else(|| Error::Solver("singular rounding transform".into()))?;
        z = &inverse * &mean;
        // Use the rotated basis, whose coordinates are near isotropic in the polytope.
        let ratio = scales.iter().fold(0.0f64, |m, v| m.max(*v)) / scales.iter().fold(f64::INFINITY, |m, v| m.min(*v));
        basis = new_basis;
        // Re-orthonormalization is not wanted: the scaled basis is the point of rounding.
        if ratio < 1.5 {
            break;
        }
    }

    // Production chains, each from the rounded centre with its own stream.
    let ids: Vec<usize> = (0..options.chains).collect();
    let thin = options.thinning.max(1) * d.max(1) / d.max(1);
    let n_rxn = model.reactions.len();
    let chains: Vec<Vec<Vec<f64>>> = crate::parallel::map_chunks(&ids, 1, |slice| {
        let mut rng = ChaCha8Rng::seed_from_u64(options.seed.wrapping_add(1_000 + slice[0] as u64));
        let mut samples = Vec::with_capacity(options.samples_per_chain);
        let total = options.burn_in + options.samples_per_chain * thin;
        let mut counter = 0usize;
        let mut collect = |v: &DVector<f64>| {
            counter += thin;
            if counter > options.burn_in && samples.len() < options.samples_per_chain {
                samples.push(lift(&pt, n_rxn, v));
            }
        };
        run(&basis, &z, &mut rng, total + thin, thin, &mut collect);
        samples
    });
    Ok(FluxSamples {
        reactions: model.reactions.iter().map(|r| r.id.clone()).collect(),
        diagnostics: diagnose(&chains, n_rxn),
        chains,
        dimension: d,
    })
}

/// Parallel artificial-centering hit-and-run.
pub fn optgp(model: &Model, options: &SamplerOptions) -> Result<FluxSamples> {
    let backend = lp::default_backend();
    let pt = build_polytope(model, options.fraction_of_optimum, backend.as_ref())?;
    let n_free = pt.free.len();
    if n_free == 0 {
        return Err(Error::Invalid("no reaction can vary; there is nothing to sample".into()));
    }
    // Warmup: the vertices that extremize each free reaction.
    let mut problem = model.to_lp();
    if let Some(f) = options.fraction_of_optimum {
        let optimum = super::fba::fba_with(model, backend.as_ref())?.objective;
        let coefs: Vec<(usize, f64)> = model.objective.clone();
        problem.add_row(f * optimum, f64::INFINITY, coefs);
    }
    problem.obj.iter_mut().for_each(|c| *c = 0.0);
    let targets: Vec<usize> = pt.free.clone();
    let parts: Vec<Result<Vec<DVector<f64>>>> =
        crate::parallel::map_chunks(&targets, crate::parallel::chunk_size(targets.len()), |slice| {
            let mut session = backend.open(&problem)?;
            let mut obj = vec![0.0; problem.n_cols()];
            let mut out = Vec::with_capacity(slice.len() * 2);
            for &j in slice {
                obj[j] = 1.0;
                for sense in [Sense::Minimize, Sense::Maximize] {
                    session.set_objective(sense, &obj)?;
                    let s = session.solve()?;
                    if s.is_optimal() {
                        out.push(DVector::from_iterator(n_free, pt.free.iter().map(|&k| s.x[k])));
                    }
                }
                obj[j] = 0.0;
            }
            Ok(out)
        });
    let mut warmup = Vec::with_capacity(targets.len() * 2);
    for part in parts {
        warmup.extend(part?);
    }
    if warmup.len() < 2 {
        return Err(Error::Solver("no warmup points were found".into()));
    }
    let center = warmup.iter().fold(DVector::zeros(n_free), |acc, p| acc + p) / warmup.len() as f64;
    let n = orthonormal_nullspace(&pt.eq);
    let d = n.ncols();
    let x0: DVector<f64> = if pt.eq.nrows() == 0 {
        DVector::zeros(n_free)
    } else {
        let gram = &pt.eq * pt.eq.transpose();
        let y = gram.lu().solve(&pt.eq_rhs).ok_or_else(|| Error::Solver("singular equality system".into()))?;
        pt.eq.transpose() * y
    };
    let project = |v: &DVector<f64>| -> DVector<f64> { &x0 + &n * (n.transpose() * (v - &x0)) };
    let ids: Vec<usize> = (0..options.chains).collect();
    let thin = options.thinning.max(1);
    let n_rxn = model.reactions.len();
    let warmup = &warmup;
    let center = &center;
    let chains: Vec<Vec<Vec<f64>>> = crate::parallel::map_chunks(&ids, 1, |slice| {
        let mut rng = ChaCha8Rng::seed_from_u64(options.seed.wrapping_add(2_000 + slice[0] as u64));
        // Every chain starts at the centre and is pushed off along one warmup direction.
        let mut v = center.clone();
        let mut samples = Vec::with_capacity(options.samples_per_chain);
        let total = options.burn_in + options.samples_per_chain * thin;
        for step in 1..=total {
            let w = &warmup[rng.random_range(0..warmup.len())];
            let dir = w - center;
            let norm = dir.norm();
            if norm > 1e-12 {
                let unit = dir / norm;
                let mut lo = f64::NEG_INFINITY;
                let mut hi = f64::INFINITY;
                for k in 0..n_free {
                    let a = unit[k];
                    if a.abs() > 1e-13 {
                        let (t1, t2) = ((pt.lb[k] - v[k]) / a, (pt.ub[k] - v[k]) / a);
                        let (tl, th) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
                        lo = lo.max(tl);
                        hi = hi.min(th);
                    }
                }
                if hi > lo {
                    let t = lo + rng.random::<f64>() * (hi - lo);
                    v.axpy(t, &unit, 1.0);
                }
            }
            if step % 100 == 0 {
                v = project(&v);
            }
            if step > options.burn_in && (step - options.burn_in) % thin == 0 {
                samples.push(lift(&pt, n_rxn, &v));
            }
        }
        samples
    });
    Ok(FluxSamples {
        reactions: model.reactions.iter().map(|r| r.id.clone()).collect(),
        diagnostics: diagnose(&chains, n_rxn),
        chains,
        dimension: d,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp::Sense;
    use crate::model::{Metabolite, Model, Reaction};

    /// `-> A`, then A leaves by one of three routes: a simplex scaled by the input.
    fn branch() -> Model {
        let mut m = Model::new("branch");
        m.metabolites.push(Metabolite { id: "A".into(), compartment: "c".into(), ..Metabolite::default() });
        let reaction = |id: &str, coef: f64, lb: f64, ub: f64| Reaction {
            id: id.into(),
            lb,
            ub,
            stoichiometry: vec![(0, coef)],
            ..Reaction::default()
        };
        m.reactions.push(reaction("in", 1.0, 1.0, 1.0));
        for k in 0..3 {
            m.reactions.push(reaction(&format!("out{k}"), -1.0, 0.0, 1.0));
        }
        m.objective = vec![(1, 1.0)];
        m.sense = Some(Sense::Maximize);
        m.reindex();
        m
    }

    /// Uniform on the 2-simplex: each coordinate has mean 1/3 and variance 1/18.
    #[test]
    fn chrr_is_uniform_on_a_simplex() {
        let options = SamplerOptions { chains: 4, samples_per_chain: 4000, thinning: 10, burn_in: 1000, seed: 3, fraction_of_optimum: None };
        let s = chrr(&branch(), &options).unwrap();
        assert_eq!(s.dimension, 2);
        for r in 1..4 {
            assert!((s.diagnostics.mean[r] - 1.0 / 3.0).abs() < 4.0 * s.diagnostics.mcse[r] + 1e-3, "{:?}", s.diagnostics.mean);
            assert!((s.diagnostics.sd[r] - (1.0f64 / 18.0).sqrt()).abs() < 0.01, "{:?}", s.diagnostics.sd);
        }
        assert!(s.diagnostics.max_split_rhat < 1.02);
    }

    #[test]
    fn optgp_keeps_samples_feasible() {
        let options = SamplerOptions { chains: 2, samples_per_chain: 500, thinning: 10, burn_in: 200, seed: 3, fraction_of_optimum: None };
        let s = optgp(&branch(), &options).unwrap();
        for chain in &s.chains {
            for v in chain {
                assert!((v[1] + v[2] + v[3] - 1.0).abs() < 1e-6);
                assert!(v[1] >= -1e-9 && v[2] >= -1e-9 && v[3] >= -1e-9);
            }
        }
    }

    #[test]
    fn diagnostics_on_known_series() {
        // Independent draws have an effective sample size close to their number.
        let mut rng = ChaCha8Rng::seed_from_u64(9);
        let chains: Vec<Vec<f64>> = (0..4).map(|_| (0..2000).map(|_| rng.random::<f64>()).collect()).collect();
        let ess = effective_sample_size(&chains);
        assert!(ess > 6000.0, "{ess}");
        assert!((split_rhat(&chains) - 1.0).abs() < 0.01);
        // Chains that disagree are flagged.
        let shifted: Vec<Vec<f64>> = chains.iter().enumerate().map(|(k, c)| c.iter().map(|v| v + k as f64).collect()).collect();
        assert!(split_rhat(&shifted) > 1.5);
    }
}
