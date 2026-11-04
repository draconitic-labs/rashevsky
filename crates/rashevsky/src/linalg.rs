//! Dense linear algebra helpers built on Gaussian elimination with full pivoting.

use nalgebra::DMatrix;

/// Rank-revealing reduced row echelon decomposition.
#[derive(Clone, Debug)]
pub struct Rref {
    /// Numerical rank.
    pub rank: usize,
    /// Pivot column of each pivot row, in row order.
    pub pivot_cols: Vec<usize>,
    /// Row indices of the original matrix chosen as independent, in pivot order.
    pub pivot_rows: Vec<usize>,
    /// The reduced matrix, with pivot rows first.
    pub reduced: DMatrix<f64>,
}

/// Reduces `a` with full pivoting. Entries below `tol` times the largest entry count as zero.
pub fn rref(a: &DMatrix<f64>, tol: f64) -> Rref {
    let (m, n) = a.shape();
    let mut r = a.clone();
    let scale = a.iter().fold(0.0f64, |acc, v| acc.max(v.abs())).max(1e-300);
    let threshold = tol * scale;
    let mut row_of: Vec<usize> = (0..m).collect();
    let mut col_free: Vec<bool> = vec![true; n];
    let mut pivot_cols = Vec::new();
    let mut rank = 0;
    while rank < m {
        let mut best = threshold;
        let mut pick = None;
        for i in rank..m {
            for j in 0..n {
                if col_free[j] && r[(i, j)].abs() > best {
                    best = r[(i, j)].abs();
                    pick = Some((i, j));
                }
            }
        }
        let Some((pi, pj)) = pick else { break };
        r.swap_rows(rank, pi);
        row_of.swap(rank, pi);
        let p = r[(rank, pj)];
        for j in 0..n {
            r[(rank, j)] /= p;
        }
        for i in 0..m {
            if i != rank {
                let f = r[(i, pj)];
                if f != 0.0 {
                    for j in 0..n {
                        let v = r[(rank, j)];
                        r[(i, j)] -= f * v;
                    }
                }
            }
        }
        col_free[pj] = false;
        pivot_cols.push(pj);
        rank += 1;
    }
    Rref { rank, pivot_cols, pivot_rows: row_of[..rank].to_vec(), reduced: r }
}

/// A basis of the right null space of `a`, one basis vector per column.
pub fn nullspace(a: &DMatrix<f64>, tol: f64) -> DMatrix<f64> {
    let n = a.ncols();
    let e = rref(a, tol);
    let mut is_pivot = vec![false; n];
    for &j in &e.pivot_cols {
        is_pivot[j] = true;
    }
    let free: Vec<usize> = (0..n).filter(|&j| !is_pivot[j]).collect();
    let mut basis = DMatrix::zeros(n, free.len());
    for (k, &f) in free.iter().enumerate() {
        basis[(f, k)] = 1.0;
        for (row, &p) in e.pivot_cols.iter().enumerate() {
            basis[(p, k)] = -e.reduced[(row, f)];
        }
    }
    basis
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nullspace_annihilates_and_has_expected_dimension() {
        let a = DMatrix::from_row_slice(3, 5, &[
            1.0, -1.0, 0.0, 0.0, 2.0,
            0.0, 1.0, -1.0, 0.0, -1.0,
            1.0, 0.0, -1.0, 0.0, 1.0, // sum of the first two rows
        ]);
        let e = rref(&a, 1e-10);
        assert_eq!(e.rank, 2);
        let n = nullspace(&a, 1e-10);
        assert_eq!(n.ncols(), 3);
        assert!((&a * &n).norm() < 1e-12);
    }
}
