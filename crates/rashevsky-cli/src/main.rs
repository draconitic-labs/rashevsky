//! `rashevsky`: batch command line over the metabolic modelling engine. Results are JSON.

use clap::{Args, Parser, Subcommand};
use nalgebra as nalgebra_matrix;
use rashevsky::cb;
use rashevsky::kinetic::{KineticModel, SimOptions, SteadyOptions, SteadyState, simulate, steady_state};
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
    /// Project a pFBA flux solution onto the Biochemistry site's map (JSON, optional SVG).
    Overlay {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        /// The site's scene file.
        #[arg(long)]
        scene: PathBuf,
        /// Rhea identifier index built by scripts/build_scene_index.py.
        #[arg(long)]
        rhea: Option<PathBuf>,
        /// Also draw the overlay to this SVG file.
        #[arg(long)]
        svg: Option<PathBuf>,
    },
    /// Write an Escher map with a force-directed layout and its reaction data from pFBA.
    Escher {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        /// Map file to write.
        #[arg(long)]
        map: PathBuf,
        /// Reaction data file to write.
        #[arg(long)]
        data: PathBuf,
        /// Only this subsystem, when the model has them.
        #[arg(long)]
        subsystem: Option<String>,
        #[arg(long, default_value_t = 1)]
        seed: u64,
    },
    /// Uniform flux samples with convergence diagnostics.
    Sample {
        model: PathBuf,
        #[command(flatten)]
        edits: Edits,
        /// `chrr` or `optgp`.
        #[arg(long, default_value = "chrr")]
        method: String,
        #[arg(long, default_value_t = 4)]
        chains: usize,
        /// Samples per chain.
        #[arg(long, default_value_t = 1000)]
        samples: usize,
        #[arg(long, default_value_t = 100)]
        thinning: usize,
        #[arg(long, default_value_t = 2000)]
        burn_in: usize,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// Keep at least this fraction of the optimal objective.
        #[arg(long)]
        fraction: Option<f64>,
        /// Include the samples in the output, not only the diagnostics.
        #[arg(long)]
        full: bool,
    },
    /// Thermodynamics-based flux analysis with energies from a JSON file.
    Tfa {
        model: PathBuf,
        energies: PathBuf,
        #[command(flatten)]
        edits: Edits,
        /// Lower concentration bound in mol/L.
        #[arg(long, default_value_t = 1e-6)]
        c_min: f64,
        /// Upper concentration bound in mol/L.
        #[arg(long, default_value_t = 0.02)]
        c_max: f64,
        /// Allowed error in standard uncertainties of the energies.
        #[arg(long, default_value_t = 2.0)]
        sigmas: f64,
    },
    /// Enzyme-constrained FBA with kcat and molecular weight data from a JSON file.
    Ecfba {
        model: PathBuf,
        enzymes: PathBuf,
        #[command(flatten)]
        edits: Edits,
        /// Replace the protein pool of the file, in g/gDW.
        #[arg(long)]
        pool: Option<f64>,
    },
    /// Time course of a kinetic model.
    Simulate {
        model: PathBuf,
        /// End time.
        #[arg(long)]
        end: f64,
        /// Number of intervals between outputs.
        #[arg(long, default_value_t = 100)]
        steps: usize,
        /// `rodas4` or `bdf`.
        #[arg(long, default_value = "rodas4")]
        method: String,
        #[arg(long, default_value_t = 1e-9)]
        rtol: f64,
        #[arg(long, default_value_t = 1e-12)]
        atol: f64,
        /// Comma-separated symbols; all species and reactions when absent.
        #[arg(long)]
        symbols: Option<String>,
        /// Override a symbol's initial value as `name=value` (repeatable).
        #[arg(long = "set", allow_hyphen_values = true)]
        set: Vec<String>,
    },
    /// Steady state of a kinetic model.
    SteadyState {
        model: PathBuf,
        /// Integrate this long before Newton iteration.
        #[arg(long, default_value_t = 100.0)]
        integrate: f64,
        #[arg(long = "set", allow_hyphen_values = true)]
        set: Vec<String>,
    },
    /// Metabolic control analysis at the steady state.
    Mca {
        model: PathBuf,
        #[arg(long, default_value_t = 100.0)]
        integrate: f64,
        #[arg(long = "set", allow_hyphen_values = true)]
        set: Vec<String>,
    },
    /// Local steady-state sensitivities of fluxes and states to constant symbols.
    Sensitivity {
        model: PathBuf,
        /// Comma-separated parameter names.
        #[arg(long)]
        parameters: String,
        #[arg(long, default_value_t = 100.0)]
        integrate: f64,
        #[arg(long = "set", allow_hyphen_values = true)]
        set: Vec<String>,
    },
    /// First-order and total Sobol indices of one output over parameter ranges.
    Sobol {
        model: PathBuf,
        /// `name=low:high`, or `name=low:high:log` to sample the logarithm (repeatable).
        #[arg(long = "range", required = true)]
        ranges: Vec<String>,
        /// Species, reaction or parameter to read.
        #[arg(long)]
        output: String,
        /// Read the output at this time rather than at the steady state.
        #[arg(long)]
        time: Option<f64>,
        #[arg(long, default_value_t = 100.0)]
        integrate: f64,
        #[arg(long, default_value_t = 512)]
        samples: usize,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        #[arg(long = "set", allow_hyphen_values = true)]
        set: Vec<String>,
    },
    /// Fits parameters to measurements by Levenberg-Marquardt, with likelihood profiles.
    Fit {
        model: PathBuf,
        /// `name=start`, with a positive start (repeatable).
        #[arg(long = "start", required = true)]
        starts: Vec<String>,
        /// CSV with the header `time,symbol,value,sigma`.
        #[arg(long)]
        data: PathBuf,
        /// Profile this fitted parameter (repeatable).
        #[arg(long = "profile")]
        profiles: Vec<String>,
        /// Profile over this factor either side of the estimate.
        #[arg(long, default_value_t = 4.0)]
        factor: f64,
        /// Grid values per profile.
        #[arg(long, default_value_t = 15)]
        points: usize,
        #[arg(long = "set", allow_hyphen_values = true)]
        set: Vec<String>,
    },
    /// Ensemble of stochastic simulations.
    Ssa {
        model: PathBuf,
        #[arg(long)]
        end: f64,
        #[arg(long, default_value_t = 100)]
        steps: usize,
        #[arg(long, default_value_t = 1000)]
        runs: usize,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// Use tau-leaping with this epsilon rather than the exact direct method.
        #[arg(long)]
        tau_epsilon: Option<f64>,
        /// Comma-separated symbols; all species when absent.
        #[arg(long)]
        symbols: Option<String>,
    },
    /// Run SBML Test Suite semantic cases and report each outcome.
    SbmlSuite {
        /// Directory that holds the numbered case folders.
        root: PathBuf,
        /// Comma-separated case numbers; all when absent.
        #[arg(long)]
        cases: Option<String>,
        /// `rodas4` or `bdf`.
        #[arg(long, default_value = "rodas4")]
        method: String,
        /// Read this SBML version, such as `l2v5`, rather than the newest in each case.
        #[arg(long)]
        level: Option<String>,
        /// Write each model and read it back before running it.
        #[arg(long)]
        round_trip: bool,
    },
}

fn load_kinetic(path: &PathBuf, set: &[String]) -> Result<KineticModel> {
    let mut km = KineticModel::read_sbml(path)?;
    for item in set {
        let invalid = || rashevsky::Error::Invalid(format!("`{item}` is not name=value"));
        let (name, value) = item.split_once('=').ok_or_else(invalid)?;
        km.set_value(name, value.parse().map_err(|_| invalid())?)?;
    }
    Ok(km)
}

fn parse_range(text: &str) -> Result<rashevsky::kinetic::sensitivity::Range> {
    let invalid = || rashevsky::Error::Invalid(format!("`{text}` is not name=low:high or name=low:high:log"));
    let (name, rest) = text.split_once('=').ok_or_else(invalid)?;
    let parts: Vec<&str> = rest.split(':').collect();
    let log = match parts.get(2) {
        None => false,
        Some(&"log") => true,
        Some(_) => return Err(invalid()),
    };
    if parts.len() < 2 || parts.len() > 3 {
        return Err(invalid());
    }
    Ok(rashevsky::kinetic::sensitivity::Range {
        name: name.to_string(),
        low: parts[0].parse().map_err(|_| invalid())?,
        high: parts[1].parse().map_err(|_| invalid())?,
        log,
    })
}

fn read_observations(path: &PathBuf) -> Result<Vec<rashevsky::kinetic::fit::Observation>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| rashevsky::Error::Invalid(format!("cannot read {}: {e}", path.display())))?;
    let mut lines = text.lines().enumerate().filter(|(_, l)| !l.trim().is_empty());
    let header: Vec<&str> = lines.next().map(|(_, l)| l.split(',').map(str::trim).collect()).unwrap_or_default();
    if header != ["time", "symbol", "value", "sigma"] {
        return Err(rashevsky::Error::Invalid(format!(
            "{}: the header must be time,symbol,value,sigma",
            path.display()
        )));
    }
    lines
        .map(|(number, line)| {
            let invalid = || rashevsky::Error::Invalid(format!("{}: line {} is malformed", path.display(), number + 1));
            let cells: Vec<&str> = line.split(',').map(str::trim).collect();
            let [time, symbol, value, sigma] = cells[..] else {
                return Err(invalid());
            };
            let sigma: f64 = sigma.parse().map_err(|_| invalid())?;
            if sigma <= 0.0 {
                return Err(invalid());
            }
            Ok(rashevsky::kinetic::fit::Observation {
                time: time.parse().map_err(|_| invalid())?,
                symbol: symbol.to_string(),
                value: value.parse().map_err(|_| invalid())?,
                sigma,
            })
        })
        .collect()
}

fn chosen_symbols(km: &KineticModel, symbols: Option<String>) -> Vec<String> {
    match symbols {
        Some(list) => list.split(',').map(|s| s.trim().to_string()).collect(),
        None => km.species.iter().map(|s| s.id.clone()).chain(km.reactions.iter().map(|r| r.id.clone())).collect(),
    }
}

fn steady_json(km: &KineticModel, state: &SteadyState) -> serde_json::Value {
    let species: BTreeMap<String, f64> = km.species.iter().map(|s| (s.id.clone(), state.values[s.view])).collect();
    let fluxes: BTreeMap<String, f64> =
        km.reactions.iter().zip(&state.fluxes).map(|(r, v)| (r.id.clone(), *v)).collect();
    serde_json::json!({
        "species": species, "fluxes": fluxes, "residual": state.residual, "newton_iterations": state.iterations,
        "time": state.time, "largest_eigenvalue": state.largest_eigenvalue, "stable": state.stable,
    })
}

fn parse_method(text: &str) -> Result<rashevsky::kinetic::Method> {
    match text {
        "rodas4" => Ok(rashevsky::kinetic::Method::Rodas4),
        "bdf" => Ok(rashevsky::kinetic::Method::Bdf),
        other => Err(rashevsky::Error::Invalid(format!("unknown method `{other}`"))),
    }
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
    } else if edits.minimize {
        model.sense = Some(Sense::Minimize);
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

/// A solution a further analysis builds on must be optimal.
fn optimal(s: cb::FluxSolution) -> Result<cb::FluxSolution> {
    if s.status == lp::Status::Optimal {
        Ok(s)
    } else {
        Err(rashevsky::Error::Solver(format!("the reference flux solve ended with {:?}", s.status)))
    }
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
    let shadow_prices = s.shadow_prices.iter().zip(&model.metabolites).map(|(v, m)| (m.id.clone(), *v)).collect();
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
                if m.is_boundary_reaction(j) {
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
                rashevsky::sbml::write::write_string(&doc)?
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
                Some(list) => Some(list.split(',').map(|id| m.reaction_index(id.trim())).collect::<Result<Vec<_>>>()?),
            };
            let ranges = cb::fva(&m, &cb::FvaOptions { fraction_of_optimum: fraction, reactions: targets })?;
            // JSON has no infinity, so an open side is reported by name as well as by a null.
            let rows: Vec<serde_json::Value> = ranges
                .iter()
                .map(|r| {
                    serde_json::json!({
                        "reaction": r.reaction, "id": r.id, "minimum": r.minimum, "maximum": r.maximum,
                        "unbounded_below": r.minimum == f64::NEG_INFINITY, "unbounded_above": r.maximum == f64::INFINITY,
                    })
                })
                .collect();
            emit(&out, &rows)
        }
        Command::Loopless { model, edits, big_m } => {
            let m = load(&model, &edits)?;
            let s = cb::loopless_fba(&m, big_m)?;
            emit(&out, &fluxes(&m, s))
        }
        Command::CycleFree { model, edits } => {
            let m = load(&model, &edits)?;
            let base = optimal(cb::fba(&m)?)?;
            let free = cb::cycle_free(&m, &base.fluxes)?;
            let s = cb::FluxSolution { fluxes: free, shadow_prices: Vec::new(), reduced_costs: Vec::new(), ..base };
            emit(&out, &fluxes(&m, s))
        }
        Command::Moma { model, edits, linear } => {
            let (reference, perturbed) = load_reference_and_perturbed(&model, &edits)?;
            let wt = optimal(cb::pfba(&reference, 1.0)?)?;
            let s = if linear { cb::linear_moma(&perturbed, &wt.fluxes)? } else { cb::moma(&perturbed, &wt.fluxes)? };
            emit(&out, &fluxes(&perturbed, s))
        }
        Command::Room { model, edits, delta, epsilon, relax } => {
            let (reference, perturbed) = load_reference_and_perturbed(&model, &edits)?;
            let wt = optimal(cb::pfba(&reference, 1.0)?)?;
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
        Command::Overlay { model, edits, scene, rhea, svg } => {
            let m = load(&model, &edits)?;
            let scene: rashevsky::io::overlay::Scene = serde_json::from_str(&std::fs::read_to_string(scene)?)?;
            let rhea_events: BTreeMap<String, Vec<String>> = match rhea {
                Some(path) => serde_json::from_str(&std::fs::read_to_string(path)?)?,
                None => BTreeMap::new(),
            };
            let s = optimal(cb::pfba(&m, 1.0)?)?;
            let objective = m.objective.iter().map(|&(j, c)| c * s.fluxes[j]).sum::<f64>();
            let overlay = rashevsky::io::overlay::build_overlay(&m, &s.fluxes, Some(objective), &scene, &rhea_events)?;
            if let Some(path) = svg {
                std::fs::write(path, rashevsky::io::overlay::render_svg(&scene, &overlay))?;
            }
            emit(&out, &overlay)
        }
        Command::Escher { model, edits, map, data, subsystem, seed } => {
            let m = load(&model, &edits)?;
            let reactions: Vec<usize> = (0..m.reactions.len())
                .filter(|&j| !m.reactions[j].is_biomass())
                .filter(|&j| subsystem.as_ref().is_none_or(|sub| m.reactions[j].subsystem.as_ref() == Some(sub)))
                .collect();
            if reactions.is_empty() {
                return Err(rashevsky::Error::Invalid(format!(
                    "no reactions to draw in {}",
                    subsystem.as_deref().unwrap_or(&m.id)
                )));
            }
            let layout = rashevsky::io::layout::force_layout(&m, &reactions, seed, 400);
            let name = subsystem.clone().unwrap_or_else(|| m.id.clone());
            let map_doc = rashevsky::io::escher::map_json(&m, &layout, &name);
            std::fs::write(map, serde_json::to_string(&map_doc)?)?;
            let s = optimal(cb::pfba(&m, 1.0)?)?;
            std::fs::write(data, serde_json::to_string_pretty(&rashevsky::io::escher::reaction_data(&m, &s.fluxes))?)?;
            emit(
                &out,
                &serde_json::json!({"nodes": layout.metabolites.len() + layout.reactions.len(), "reactions": layout.reactions.len(), "objective": m.objective.iter().map(|&(j, c)| c * s.fluxes[j]).sum::<f64>()}),
            )
        }
        Command::Sample { model, edits, method, chains, samples, thinning, burn_in, seed, fraction, full } => {
            let m = load(&model, &edits)?;
            let options = cb::SamplerOptions {
                chains,
                samples_per_chain: samples,
                thinning,
                burn_in,
                seed,
                fraction_of_optimum: fraction,
            };
            let started = std::time::Instant::now();
            let result = match method.as_str() {
                "chrr" => cb::chrr(&m, &options)?,
                "optgp" => cb::optgp(&m, &options)?,
                other => return Err(rashevsky::Error::Invalid(format!("unknown sampler `{other}`"))),
            };
            let seconds = started.elapsed().as_secs_f64();
            if full {
                emit(&out, &result)
            } else {
                emit(
                    &out,
                    &serde_json::json!({
                        "reactions": result.reactions,
                        "dimension": result.dimension,
                        "diagnostics": result.diagnostics,
                        "seconds": seconds,
                    }),
                )
            }
        }
        Command::Tfa { model, energies, edits, c_min, c_max, sigmas } => {
            let m = load(&model, &edits)?;
            let data: cb::ThermoData = serde_json::from_str(&std::fs::read_to_string(energies)?)?;
            let options = cb::TfaOptions { c_min, c_max, sigmas, ..cb::TfaOptions::default() };
            let s = cb::tfa(&m, &data, &options)?;
            emit(
                &out,
                &serde_json::json!({
                    "fluxes": fluxes(&m, s.flux),
                    "concentrations": s.concentrations,
                    "gibbs": s.gibbs,
                }),
            )
        }
        Command::Ecfba { model, enzymes, edits, pool } => {
            let m = load(&model, &edits)?;
            let mut data: cb::EnzymeData = serde_json::from_str(&std::fs::read_to_string(enzymes)?)?;
            if let Some(p) = pool {
                data.protein_pool = p;
            }
            let s = cb::enzyme_constrained_fba(&m, &data)?;
            emit(
                &out,
                &serde_json::json!({
                    "fluxes": fluxes(&m, s.flux),
                    "usage": s.usage,
                    "pool_used": s.pool_used,
                }),
            )
        }
        Command::Simulate { model, end, steps, method, rtol, atol, symbols, set } => {
            let km = load_kinetic(&model, &set)?;
            let times: Vec<f64> = (0..=steps).map(|k| end * k as f64 / steps as f64).collect();
            let options = SimOptions { method: parse_method(&method)?, rtol, atol, ..SimOptions::default() };
            let course = simulate(&km, &times, &options)?;
            let chosen = chosen_symbols(&km, symbols);
            let mut series = BTreeMap::new();
            for name in &chosen {
                series.insert(name.clone(), course.series(&km, name)?);
            }
            emit(
                &out,
                &serde_json::json!({
                    "time": course.time,
                    "series": series,
                    "stats": {"steps": course.stats.steps, "rejected": course.stats.rejected, "rhs_evals": course.stats.rhs_evals, "jacobians": course.stats.jacobians},
                }),
            )
        }
        Command::SteadyState { model, integrate, set } => {
            let km = load_kinetic(&model, &set)?;
            let state = steady_state(&km, &SteadyOptions { integrate_first: integrate, ..SteadyOptions::default() })?;
            emit(&out, &steady_json(&km, &state))
        }
        Command::Mca { model, integrate, set } => {
            let km = load_kinetic(&model, &set)?;
            let state = steady_state(&km, &SteadyOptions { integrate_first: integrate, ..SteadyOptions::default() })?;
            let c = rashevsky::kinetic::control_analysis(&km, &state)?;
            let rows = |m: &nalgebra_matrix::DMatrix<f64>| -> Vec<Vec<f64>> {
                (0..m.nrows()).map(|i| (0..m.ncols()).map(|j| m[(i, j)]).collect()).collect()
            };
            emit(
                &out,
                &serde_json::json!({
                    "steady_state": steady_json(&km, &state),
                    "reactions": c.reactions,
                    "states": c.states,
                    "flux_control": rows(&c.scaled_flux_control),
                    "concentration_control": rows(&c.scaled_concentration_control),
                    "elasticities": rows(&c.scaled_elasticities),
                    "theorems": c.theorems,
                    "conserved_moieties": c.states.len() - c.moieties.independent.len(),
                }),
            )
        }
        Command::Sensitivity { model, parameters, integrate, set } => {
            let km = load_kinetic(&model, &set)?;
            let state = steady_state(&km, &SteadyOptions { integrate_first: integrate, ..SteadyOptions::default() })?;
            let control = rashevsky::kinetic::control_analysis(&km, &state)?;
            let names: Vec<String> = parameters.split(',').map(|p| p.trim().to_string()).collect();
            let local = rashevsky::kinetic::sensitivity::local_steady(&km, &state, &control, &names)?;
            let rows = |m: &nalgebra_matrix::DMatrix<f64>| -> Vec<Vec<f64>> {
                (0..m.nrows()).map(|i| (0..m.ncols()).map(|j| m[(i, j)]).collect()).collect()
            };
            emit(
                &out,
                &serde_json::json!({
                    "parameters": local.parameters,
                    "reactions": control.reactions,
                    "states": control.states,
                    "flux": rows(&local.flux),
                    "state": rows(&local.state),
                }),
            )
        }
        Command::Sobol { model, ranges, output, time, integrate, samples, seed, set } => {
            use rashevsky::kinetic::sensitivity::{Output, sobol_model};
            let km = load_kinetic(&model, &set)?;
            let ranges = ranges.iter().map(|r| parse_range(r)).collect::<Result<Vec<_>>>()?;
            let reading = match time {
                Some(time) => Output::AtTime { symbol: output, time },
                None => Output::Steady { symbol: output, integrate },
            };
            emit(&out, &sobol_model(&km, &ranges, samples, seed, &reading)?)
        }
        Command::Fit { model, starts, data, profiles, factor, points, set } => {
            use rashevsky::kinetic::fit::{FitParameter, fit, profile};
            let km = load_kinetic(&model, &set)?;
            let parameters = starts
                .iter()
                .map(|item| {
                    let invalid = || rashevsky::Error::Invalid(format!("`{item}` is not name=start"));
                    let (name, start) = item.split_once('=').ok_or_else(invalid)?;
                    Ok(FitParameter { name: name.to_string(), start: start.parse().map_err(|_| invalid())? })
                })
                .collect::<Result<Vec<_>>>()?;
            let observations = read_observations(&data)?;
            let options = rashevsky::kinetic::SimOptions::default();
            let fitted = fit(&km, &parameters, &observations, &options)?;
            let curves = profiles
                .iter()
                .map(|name| profile(&km, &parameters, &observations, &fitted, name, factor, points, &options))
                .collect::<Result<Vec<_>>>()?;
            emit(&out, &serde_json::json!({"fit": fitted, "profiles": curves}))
        }
        Command::Ssa { model, end, steps, runs, seed, tau_epsilon, symbols } => {
            let km = KineticModel::read_sbml(&model)?;
            let times: Vec<f64> = (0..=steps).map(|k| end * k as f64 / steps as f64).collect();
            let algorithm = match tau_epsilon {
                Some(epsilon) => rashevsky::kinetic::ssa::Algorithm::TauLeap { epsilon },
                None => rashevsky::kinetic::ssa::Algorithm::Direct,
            };
            let e = rashevsky::kinetic::ssa::ensemble(&km, &times, runs, algorithm, seed)?;
            let chosen = chosen_symbols(&km, symbols);
            let mut mean = BTreeMap::new();
            let mut sd = BTreeMap::new();
            for name in &chosen {
                let slot = km.slot(name)?;
                mean.insert(name.clone(), e.mean.iter().map(|v| v[slot]).collect::<Vec<f64>>());
                sd.insert(name.clone(), e.sd.iter().map(|v| v[slot]).collect::<Vec<f64>>());
            }
            emit(&out, &serde_json::json!({"time": e.time, "runs": e.runs, "mean": mean, "sd": sd}))
        }
        Command::SbmlSuite { root, cases: list, method, level, round_trip } => {
            let filter: Option<Vec<String>> = list.map(|l| l.split(',').map(|c| format!("{:0>5}", c.trim())).collect());
            let options = rashevsky::suite::SuiteOptions { method: parse_method(&method)?, level, round_trip };
            emit(&out, &rashevsky::suite::run_suite(&root, &options, filter.as_deref())?)
        }
    }
}

fn main() {
    // Deeply nested SBML math recurses in the evaluator, so workers get a large stack.
    let _ = rayon::ThreadPoolBuilder::new().stack_size(64 << 20).build_global();
    if let Err(e) = run(Cli::parse()) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}
