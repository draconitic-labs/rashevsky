//! `rashevsky`: batch command line over the metabolic modelling engine. Results are JSON.

mod cases;

use clap::{Args, Parser, Subcommand};
use rashevsky::cb;
use rashevsky::lp::{self, Backend, Sense};
use rashevsky::{Model, Result};
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "rashevsky", version, about = "Genome-scale constraint-based and kinetic metabolic modelling")]
struct Cli {
    /// Write the result to this file instead of standard output.
    #[arg(long, global = true)]
    out: Option<PathBuf>,
    #[command(subcommand)]
    command: Command,
}

/// Edits applied to a model before an analysis.
#[derive(Args, Clone, Default)]
struct Edits {
    /// Delete a gene (repeatable).
    #[arg(long = "ko-gene")]
    ko_gene: Vec<String>,
    /// Delete a reaction (repeatable).
    #[arg(long = "ko-reaction")]
    ko_reaction: Vec<String>,
    /// Set a reaction's bounds as `id=lb:ub` (repeatable).
    #[arg(long = "bound", allow_hyphen_values = true)]
    bound: Vec<String>,
    /// Replace the objective with this reaction.
    #[arg(long)]
    objective: Option<String>,
    /// Minimize the objective rather than maximize it.
    #[arg(long)]
    minimize: bool,
}

#[derive(Subcommand)]
enum Command {
    /// Size, compartments and objective of a model.
    Info { model: PathBuf },
    /// Mass and charge balance of every non-boundary reaction.
    Balance { model: PathBuf },
    /// Rewrite an SBML file as SBML Level 3 with FBC 2 and groups, or as JSON (by extension).
    Convert { input: PathBuf, output: PathBuf },
    /// Flux balance analysis.
    Fba {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        /// Solver: `highs` or `simplex`.
        #[arg(long, default_value = "highs")]
        backend: String,
    },
    /// Parsimonious FBA.
    Pfba {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        #[arg(long, default_value_t = 1.0)]
        fraction: f64,
    },
    /// Flux variability analysis.
    Fva {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        #[arg(long, default_value_t = 1.0)]
        fraction: f64,
        /// Comma-separated reactions; all when absent.
        #[arg(long)]
        reactions: Option<String>,
    },
    /// Mixed-integer loopless FBA.
    Loopless {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        #[arg(long, default_value_t = 1000.0)]
        big_m: f64,
    },
    /// Remove loops from the FBA solution (CycleFreeFlux).
    CycleFree {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
    },
    /// MOMA for a perturbation, with the unperturbed pFBA solution as reference.
    Moma {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        /// Minimize the L1 instead of the L2 distance.
        #[arg(long)]
        linear: bool,
    },
    /// ROOM for a perturbation, with the unperturbed pFBA solution as reference.
    Room {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        #[arg(long, default_value_t = 0.03)]
        delta: f64,
        #[arg(long, default_value_t = 0.001)]
        epsilon: f64,
        #[arg(long)]
        relax: bool,
    },
    /// Single or double gene or reaction deletions.
    Deletion {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        /// `gene` or `reaction`.
        #[arg(long, default_value = "gene")]
        kind: String,
        /// 1 or 2.
        #[arg(long, default_value_t = 1)]
        order: u8,
    },
    /// Synthetic lethal gene pairs.
    SyntheticLethal {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        /// Growth at or below this value counts as lethal.
        #[arg(long, default_value_t = 1e-3)]
        threshold: f64,
    },
    /// Production envelope of a target flux against a control flux.
    Envelope {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        #[arg(long)]
        control: String,
        #[arg(long)]
        target: String,
        #[arg(long, default_value_t = 20)]
        points: usize,
    },
    /// Minimal medium for a required growth rate.
    MinimalMedium {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        #[arg(long)]
        growth: f64,
        /// Minimize the number of components before total uptake.
        #[arg(long)]
        components: bool,
    },
    /// Dynamic FBA from a JSON settings file.
    Dfba {
        model: PathBuf,
        settings: PathBuf,
        #[command(flatten)]
        edits: Edits,
    },
    /// Evaluate a metabolic task list (RAVEN and Human-GEM format).
    Tasks {
        model: PathBuf,
        tasks: PathBuf,
        /// Keep the model's own boundary reactions open.
        #[arg(long)]
        open_boundaries: bool,
    },
    /// Run a study described by a case file (see docs/CASES.md).
    Run { case: PathBuf },
}

fn parse_bound(text: &str) -> Result<(String, f64, f64)> {
    let invalid = || rashevsky::Error::Invalid(format!("bound `{text}` is not id=lb:ub"));
    let (id, range) = text.split_once('=').ok_or_else(invalid)?;
    let (lb, ub) = range.split_once(':').ok_or_else(invalid)?;
    let parse = |s: &str| -> Result<f64> {
        match s {
            "inf" | "INF" => Ok(f64::INFINITY),
            "-inf" | "-INF" => Ok(f64::NEG_INFINITY),
            _ => s.parse().map_err(|_| invalid()),
        }
    };
    Ok((id.to_string(), parse(lb)?, parse(ub)?))
}

/// Applies bounds, objective and deletions.
fn apply(model: &mut Model, edits: &Edits) -> Result<()> {
    for b in &edits.bound {
        let (id, lb, ub) = parse_bound(b)?;
        model.set_bounds(&id, lb, ub)?;
    }
    if let Some(obj) = &edits.objective {
        model.set_objective(obj, if edits.minimize { Sense::Minimize } else { Sense::Maximize })?;
    }
    for id in &edits.ko_reaction {
        model.set_bounds(id, 0.0, 0.0)?;
    }
    if !edits.ko_gene.is_empty() {
        let mut deleted = vec![false; model.genes.len()];
        for g in &edits.ko_gene {
            deleted[model.gene_index(g)?] = true;
        }
        for j in model.reactions_disabled_by(&deleted) {
            model.reactions[j].lb = 0.0;
            model.reactions[j].ub = 0.0;
        }
    }
    Ok(())
}

fn load(path: &PathBuf, edits: &Edits) -> Result<Model> {
    let mut model = Model::read_sbml(path)?;
    apply(&mut model, edits)?;
    Ok(model)
}

/// Editing bounds only, used for the reference state of MOMA and ROOM.
fn load_reference_and_perturbed(path: &PathBuf, edits: &Edits) -> Result<(Model, Model)> {
    let mut reference = Model::read_sbml(path)?;
    let base = Edits { ko_gene: Vec::new(), ko_reaction: Vec::new(), ..edits.clone() };
    apply(&mut reference, &base)?;
    let mut perturbed = reference.clone();
    let knockouts = Edits { bound: Vec::new(), objective: None, ..edits.clone() };
    apply(&mut perturbed, &knockouts)?;
    Ok((reference, perturbed))
}

#[derive(Serialize)]
struct Fluxes {
    status: lp::Status,
    objective: f64,
    fluxes: BTreeMap<String, f64>,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    shadow_prices: BTreeMap<String, f64>,
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    reduced_costs: BTreeMap<String, f64>,
}

fn fluxes(model: &Model, s: cb::FluxSolution) -> Fluxes {
    let by_rxn = |values: &[f64]| -> BTreeMap<String, f64> {
        values.iter().zip(&model.reactions).map(|(v, r)| (r.id.clone(), *v)).collect()
    };
    let shadow_prices = s
        .shadow_prices
        .iter()
        .zip(&model.metabolites)
        .map(|(v, m)| (m.id.clone(), *v))
        .collect();
    Fluxes {
        status: s.status,
        objective: s.objective,
        fluxes: by_rxn(&s.fluxes),
        shadow_prices,
        reduced_costs: by_rxn(&s.reduced_costs),
    }
}

fn emit(out: &Option<PathBuf>, value: &impl Serialize) -> Result<()> {
    let text = serde_json::to_string_pretty(value)?;
    match out {
        Some(path) => std::fs::write(path, text + "\n")?,
        None => println!("{text}"),
    }
    Ok(())
}

fn run(cli: Cli) -> Result<()> {
    let out = cli.out.clone();
    match cli.command {
        Command::Info { model } => emit(&out, &Model::read_sbml(model)?.summary()),
        Command::Balance { model } => {
            let m = Model::read_sbml(model)?;
            let mut report = BTreeMap::new();
            let mut unknown = Vec::new();
            let mut balanced = 0usize;
            for (j, r) in m.reactions.iter().enumerate() {
                if r.is_boundary() {
                    continue;
                }
                match m.imbalance(j) {
                    None => unknown.push(r.id.clone()),
                    Some(map) if map.is_empty() => balanced += 1,
                    Some(map) => {
                        report.insert(r.id.clone(), map);
                    }
                }
            }
            emit(&out, &serde_json::json!({"balanced": balanced, "formula_missing": unknown, "imbalanced": report}))
        }
        Command::Convert { input, output } => {
            let doc = rashevsky::sbml::read::read_file(&input)?;
            let is_json = output.extension().is_some_and(|e| e == "json");
            let text = if is_json {
                serde_json::to_string_pretty(&Model::from_document(&doc)?)?
            } else {
                rashevsky::sbml::write::write_string(&doc)
            };
            std::fs::write(&output, text)?;
            Ok(())
        }
        Command::Fba { model, edits, backend } => {
            let m = load(&model, &edits)?;
            let chosen: Box<dyn Backend> = match backend.as_str() {
                "simplex" => Box::new(lp::simplex::Simplex::default()),
                "highs" => lp::default_backend(),
                other => return Err(rashevsky::Error::Invalid(format!("unknown backend `{other}`"))),
            };
            let s = cb::fba_with(&m, chosen.as_ref())?;
            emit(&out, &fluxes(&m, s))
        }
        Command::Pfba { model, edits, fraction } => {
            let m = load(&model, &edits)?;
            let s = cb::pfba(&m, fraction)?;
            emit(&out, &fluxes(&m, s))
        }
        Command::Fva { model, edits, fraction, reactions } => {
            let m = load(&model, &edits)?;
            let targets = match reactions {
                None => None,
                Some(list) => Some(
                    list.split(',')
                        .map(|id| m.reaction_index(id.trim()))
                        .collect::<Result<Vec<_>>>()?,
                ),
            };
            let ranges = cb::fva(&m, &cb::FvaOptions { fraction_of_optimum: fraction, reactions: targets })?;
            emit(&out, &ranges)
        }
        Command::Loopless { model, edits, big_m } => {
            let m = load(&model, &edits)?;
            let s = cb::loopless_fba(&m, big_m)?;
            emit(&out, &fluxes(&m, s))
        }
        Command::CycleFree { model, edits } => {
            let m = load(&model, &edits)?;
            let base = cb::fba(&m)?;
            let free = cb::cycle_free(&m, &base.fluxes)?;
            let s = cb::FluxSolution { fluxes: free, shadow_prices: Vec::new(), reduced_costs: Vec::new(), ..base };
            emit(&out, &fluxes(&m, s))
        }
        Command::Moma { model, edits, linear } => {
            let (reference, perturbed) = load_reference_and_perturbed(&model, &edits)?;
            let wt = cb::pfba(&reference, 1.0)?;
            let s = if linear {
                cb::linear_moma(&perturbed, &wt.fluxes)?
            } else {
                cb::moma(&perturbed, &wt.fluxes)?
            };
            emit(&out, &fluxes(&perturbed, s))
        }
        Command::Room { model, edits, delta, epsilon, relax } => {
            let (reference, perturbed) = load_reference_and_perturbed(&model, &edits)?;
            let wt = cb::pfba(&reference, 1.0)?;
            let s = cb::room(&perturbed, &wt.fluxes, cb::RoomOptions { delta, epsilon, relax })?;
            emit(&out, &fluxes(&perturbed, s))
        }
        Command::Deletion { model, edits, kind, order } => {
            let m = load(&model, &edits)?;
            let result = match (kind.as_str(), order) {
                ("gene", 1) => cb::single_gene_deletion(&m, None)?,
                ("gene", 2) => cb::double_gene_deletion(&m, None)?,
                ("reaction", 1) => cb::single_reaction_deletion(&m, None)?,
                ("reaction", 2) => cb::double_reaction_deletion(&m, None)?,
                _ => return Err(rashevsky::Error::Invalid("kind is gene or reaction, order 1 or 2".into())),
            };
            emit(&out, &result)
        }
        Command::SyntheticLethal { model, edits, threshold } => {
            let m = load(&model, &edits)?;
            emit(&out, &cb::synthetic_lethal_genes(&m, threshold)?)
        }
        Command::Envelope { model, edits, control, target, points } => {
            let m = load(&model, &edits)?;
            emit(&out, &cb::production_envelope(&m, &control, &target, points, None)?)
        }
        Command::MinimalMedium { model, edits, growth, components } => {
            let m = load(&model, &edits)?;
            emit(&out, &cb::minimal_medium(&m, growth, components)?)
        }
        Command::Dfba { model, settings, edits } => {
            let m = load(&model, &edits)?;
            let options: cb::DfbaOptions = serde_json::from_str(&std::fs::read_to_string(settings)?)?;
            emit(&out, &cb::dfba(&m, &options)?)
        }
        Command::Tasks { model, tasks, open_boundaries } => {
            let m = Model::read_sbml(model)?;
            let list = rashevsky::tasks::parse_task_list(&std::fs::read_to_string(tasks)?)?;
            emit(&out, &rashevsky::tasks::check_tasks(&m, &list, !open_boundaries)?)
        }
        Command::Run { case } => cases::run(&case, &out),
    }
}

fn main() {
    if let Err(e) = run(Cli::parse()) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}
