//! Single and double gene or reaction deletions, and synthetic lethal pairs.
//!
//! A gene deletion disables each reaction whose rule evaluates to false (Thiele and Palsson 2010).
//! Workers keep one solver session and restart from the previous basis after each change of bounds.

use crate::error::Result;
use crate::lp::{self, Backend, Session, Status};
use crate::model::Model;
use serde::{Deserialize, Serialize};

/// Objective value after a deletion.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Deletion {
    /// Identifiers of the deleted genes or reactions.
    pub ids: Vec<String>,
    /// Solver status.
    pub status: Status,
    /// Objective value; 0 when the deletion makes the model infeasible.
    pub growth: f64,
}

fn solve_with_disabled(model: &Model, session: &mut dyn Session, disabled: &[usize]) -> Result<(Status, f64)> {
    for &j in disabled {
        session.set_col_bounds(j, 0.0, 0.0)?;
    }
    let s = session.solve()?;
    let result = match s.status {
        Status::Optimal => (s.status, s.objective),
        other => (other, 0.0),
    };
    for &j in disabled {
        session.set_col_bounds(j, model.reactions[j].lb, model.reactions[j].ub)?;
    }
    Ok(result)
}

fn run_sets(model: &Model, sets: &[(Vec<String>, Vec<usize>)], backend: &dyn Backend) -> Result<Vec<Deletion>> {
    let problem = model.to_lp();
    let chunk = crate::parallel::chunk_size(sets.len());
    let parts: Vec<Result<Vec<Deletion>>> = crate::parallel::map_chunks(sets, chunk, |slice| {
        let mut session = backend.open(&problem)?;
        slice
            .iter()
            .map(|(ids, disabled)| {
                let (status, growth) = solve_with_disabled(model, session.as_mut(), disabled)?;
                Ok(Deletion { ids: ids.clone(), status, growth })
            })
            .collect()
    });
    let mut out = Vec::with_capacity(sets.len());
    for part in parts {
        out.extend(part?);
    }
    Ok(out)
}

fn gene_sets(model: &Model, genes: &[Vec<usize>]) -> Vec<(Vec<String>, Vec<usize>)> {
    genes
        .iter()
        .map(|set| {
            let mut deleted = vec![false; model.genes.len()];
            for &g in set {
                deleted[g] = true;
            }
            let ids = set.iter().map(|&g| model.genes[g].id.clone()).collect();
            (ids, model.reactions_disabled_by(&deleted))
        })
        .collect()
}

/// Deletes each gene in turn; `genes = None` means every gene.
pub fn single_gene_deletion(model: &Model, genes: Option<&[usize]>) -> Result<Vec<Deletion>> {
    single_gene_deletion_with(model, genes, lp::default_backend().as_ref())
}

/// [`single_gene_deletion`] with a chosen backend.
pub fn single_gene_deletion_with(
    model: &Model,
    genes: Option<&[usize]>,
    backend: &dyn Backend,
) -> Result<Vec<Deletion>> {
    let sets: Vec<Vec<usize>> = match genes {
        Some(list) => list.iter().map(|&g| vec![g]).collect(),
        None => (0..model.genes.len()).map(|g| vec![g]).collect(),
    };
    run_sets(model, &gene_sets(model, &sets), backend)
}

/// Deletes each reaction in turn; `reactions = None` means every reaction.
pub fn single_reaction_deletion(model: &Model, reactions: Option<&[usize]>) -> Result<Vec<Deletion>> {
    let backend = lp::default_backend();
    let list: Vec<usize> = reactions.map_or_else(|| (0..model.reactions.len()).collect(), <[usize]>::to_vec);
    let sets: Vec<(Vec<String>, Vec<usize>)> =
        list.iter().map(|&j| (vec![model.reactions[j].id.clone()], vec![j])).collect();
    run_sets(model, &sets, backend.as_ref())
}

/// Deletes every unordered pair of the given genes.
pub fn double_gene_deletion(model: &Model, genes: Option<&[usize]>) -> Result<Vec<Deletion>> {
    let backend = lp::default_backend();
    let list: Vec<usize> = genes.map_or_else(|| (0..model.genes.len()).collect(), <[usize]>::to_vec);
    let mut sets = Vec::with_capacity(list.len() * list.len() / 2);
    for (k, &a) in list.iter().enumerate() {
        for &b in &list[k + 1..] {
            sets.push(vec![a, b]);
        }
    }
    run_sets(model, &gene_sets(model, &sets), backend.as_ref())
}

/// Deletes every unordered pair of the given reactions.
pub fn double_reaction_deletion(model: &Model, reactions: Option<&[usize]>) -> Result<Vec<Deletion>> {
    let backend = lp::default_backend();
    let list: Vec<usize> = reactions.map_or_else(|| (0..model.reactions.len()).collect(), <[usize]>::to_vec);
    let mut sets = Vec::with_capacity(list.len() * list.len() / 2);
    for (k, &a) in list.iter().enumerate() {
        for &b in &list[k + 1..] {
            sets.push((vec![model.reactions[a].id.clone(), model.reactions[b].id.clone()], vec![a, b]));
        }
    }
    run_sets(model, &sets, backend.as_ref())
}

/// Synthetic lethal gene pairs: each single deletion leaves growth above `threshold`, the pair does not.
pub fn synthetic_lethal_genes(model: &Model, threshold: f64) -> Result<Vec<Deletion>> {
    let singles = single_gene_deletion(model, None)?;
    let viable: Vec<usize> = singles.iter().enumerate().filter(|(_, d)| d.growth > threshold).map(|(g, _)| g).collect();
    let doubles = double_gene_deletion(model, Some(&viable))?;
    Ok(doubles.into_iter().filter(|d| d.growth <= threshold).collect())
}
