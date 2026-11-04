//! Commands the interface invokes.

use crate::jobs::{JobStatus, Jobs};
use crate::payload::{ModelView, computed_layout, model_view};
use rashevsky::Model;
use rashevsky::cb;
use rashevsky::kinetic::{self, KineticModel, Method, SimOptions, SteadyOptions, ssa};
use rashevsky::lp::{self, Sense, Status};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::Ordering;
use std::sync::{Arc, Mutex};
use tauri::State;

const CORE_XML: &[u8] = include_bytes!("../../../models/e_coli_core/e_coli_core.xml.gz");
const CORE_LAYOUT: &str = include_str!("../../../models/e_coli_core/layout.json");
const TEUSINK_XML: &str = include_str!("../../../models/teusink/BIOMD0000000064.xml");
const SCENE: &str = include_str!("../../../data/biochemistry/scene.json");
const RHEA: &str = include_str!("../../../data/biochemistry/rhea-events.json");
const VALIDATION: &[(&str, &str)] = &[
    ("v01_core_cobrapy", include_str!("../../../validation/results/v01_core_cobrapy.json")),
    ("v02_iml1515_keio", include_str!("../../../validation/results/v02_iml1515_keio.json")),
    ("v03_human_gem", include_str!("../../../validation/results/v03_human_gem.json")),
    ("v04_teusink_mca", include_str!("../../../validation/results/v04_teusink_mca.json")),
    ("v05_sbml_suite", include_str!("../../../validation/results/v05_sbml_suite.json")),
    ("v06_sampling", include_str!("../../../validation/results/v06_sampling.json")),
    ("v07_overlay_escher", include_str!("../../../validation/results/v07_overlay_escher.json")),
    ("v08_numerics", include_str!("../../../validation/results/v08_numerics.json")),
    ("v11_recon3d", include_str!("../../../validation/results/v11_recon3d.json")),
    ("v12_ecmodel", include_str!("../../../validation/results/v12_ecmodel.json")),
];

/// Everything the window shares with this process.
#[derive(Default)]
pub struct AppState {
    model: Mutex<Option<Arc<Model>>>,
    kinetic: Mutex<Option<Arc<KineticModel>>>,
    jobs: Jobs,
}

type Reply<T> = Result<T, String>;

fn text(e: impl std::fmt::Display) -> String {
    e.to_string()
}

fn current(state: &State<'_, AppState>) -> Reply<Arc<Model>> {
    state.model.lock().map_err(text)?.clone().ok_or_else(|| "open a model first".to_string())
}

fn current_kinetic(state: &State<'_, AppState>) -> Reply<Arc<KineticModel>> {
    state.kinetic.lock().map_err(text)?.clone().ok_or_else(|| "open a kinetic model first".to_string())
}

fn data_dir() -> PathBuf {
    let base = std::env::var_os("XDG_DATA_HOME")
        .map(PathBuf::from)
        .or_else(|| std::env::var_os("HOME").map(|h| PathBuf::from(h).join(".local/share")))
        .unwrap_or_else(std::env::temp_dir);
    base.join("rashevsky/models")
}

/// Arguments the window starts from: an optional workspace and a model file.
#[tauri::command]
pub fn startup() -> Value {
    let mut view = None;
    let mut file = None;
    for arg in std::env::args().skip(1) {
        if let Some(v) = arg.strip_prefix("--view=") {
            view = Some(v.to_string());
        } else if !arg.starts_with('-') {
            file = Some(arg);
        }
    }
    let qa = std::env::var_os("RASHEVSKY_QA").is_some();
    json!({"view": view, "file": file, "qa": qa, "models_dir": data_dir()})
}

/// Opens a model bundled with the application.
#[tauri::command]
pub async fn open_bundled(state: State<'_, AppState>, name: String) -> Reply<ModelView> {
    if name != "e_coli_core" {
        return Err(format!("no bundled model named {name}"));
    }
    let doc = rashevsky::sbml::read::read_bytes(CORE_XML).map_err(text)?;
    let model = Model::from_document(&doc).map_err(text)?;
    let view = model_view(&model, true);
    *state.model.lock().map_err(text)? = Some(Arc::new(model));
    Ok(view)
}

/// Opens an SBML file as a constraint-based model.
#[tauri::command]
pub async fn open_file(state: State<'_, AppState>, path: String) -> Reply<ModelView> {
    let path2 = path.clone();
    let (model, view) = tauri::async_runtime::spawn_blocking(move || {
        let model = Model::read_sbml(&path2).map_err(text)?;
        let view = model_view(&model, false);
        Ok::<_, String>((model, view))
    })
    .await
    .map_err(text)??;
    *state.model.lock().map_err(text)? = Some(Arc::new(model));
    Ok(view)
}

const SOURCES: &[(&str, &str, &str)] = &[
    ("iML1515", "http://bigg.ucsd.edu/static/models/iML1515.xml.gz", "iML1515.xml.gz"),
    ("iJO1366", "http://bigg.ucsd.edu/static/models/iJO1366.xml.gz", "iJO1366.xml.gz"),
    ("Recon3D", "http://bigg.ucsd.edu/static/models/Recon3D.xml.gz", "Recon3D.xml.gz"),
    (
        "Human-GEM",
        "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/Human-GEM.xml",
        "Human-GEM.xml",
    ),
];

/// Downloads a published model once into the user's data folder and returns its path.
#[tauri::command]
pub async fn download_model(name: String) -> Reply<String> {
    tauri::async_runtime::spawn_blocking(move || {
        let (_, url, file) = SOURCES.iter().find(|s| s.0 == name).ok_or_else(|| format!("unknown model {name}"))?;
        let dir = data_dir();
        std::fs::create_dir_all(&dir).map_err(text)?;
        let path = dir.join(file);
        if path.exists() {
            return Ok(path.display().to_string());
        }
        let mut reply = ureq::get(*url).call().map_err(|e| format!("could not download {name}: {e}"))?;
        let temporary = dir.join(format!("{file}.part"));
        let mut out = std::fs::File::create(&temporary).map_err(text)?;
        std::io::copy(&mut reply.body_mut().as_reader(), &mut out).map_err(text)?;
        std::fs::rename(&temporary, &path).map_err(text)?;
        Ok(path.display().to_string())
    })
    .await
    .map_err(text)?
}

/// Perturbation of the open model.
#[derive(Clone, Debug, Default, Deserialize)]
pub struct Perturbation {
    #[serde(default)]
    pub genes: Vec<String>,
    #[serde(default)]
    pub reactions: Vec<String>,
    #[serde(default)]
    pub bounds: Vec<(String, f64, f64)>,
    #[serde(default)]
    pub objective: Option<String>,
}

fn apply_medium_and_objective(model: &Model, p: &Perturbation) -> Reply<Model> {
    let mut m = model.clone();
    for (id, lb, ub) in &p.bounds {
        m.set_bounds(id, *lb, *ub).map_err(text)?;
    }
    if let Some(obj) = &p.objective {
        m.set_objective(obj, Sense::Maximize).map_err(text)?;
    }
    Ok(m)
}

fn apply_knockouts(m: &mut Model, p: &Perturbation) -> Reply<()> {
    for id in &p.reactions {
        m.set_bounds(id, 0.0, 0.0).map_err(text)?;
    }
    if !p.genes.is_empty() {
        let mut deleted = vec![false; m.genes.len()];
        for g in &p.genes {
            deleted[m.gene_index(g).map_err(text)?] = true;
        }
        for j in m.reactions_disabled_by(&deleted) {
            m.reactions[j].lb = 0.0;
            m.reactions[j].ub = 0.0;
        }
    }
    Ok(())
}

/// A solve request.
#[derive(Deserialize)]
pub struct SolveRequest {
    pub method: String,
    #[serde(flatten)]
    pub perturbation: Perturbation,
}

/// The fluxes of a solve, in the model's reaction order.
#[derive(Serialize)]
pub struct SolveReply {
    status: Status,
    objective: Option<f64>,
    growth: Option<f64>,
    fluxes: Vec<f64>,
    disabled: Vec<String>,
    milliseconds: f64,
    total_flux: f64,
}

/// Runs one flux analysis.
#[tauri::command]
pub async fn solve(state: State<'_, AppState>, request: SolveRequest) -> Reply<SolveReply> {
    let base = current(&state)?;
    tauri::async_runtime::spawn_blocking(move || {
        let started = std::time::Instant::now();
        let reference = apply_medium_and_objective(&base, &request.perturbation)?;
        let mut perturbed = reference.clone();
        apply_knockouts(&mut perturbed, &request.perturbation)?;
        let s = match request.method.as_str() {
            "fba" => cb::fba(&perturbed),
            "pfba" => cb::pfba(&perturbed, 1.0),
            "moma" => cb::pfba(&reference, 1.0).and_then(|w| cb::moma(&perturbed, &w.fluxes)),
            "lmoma" => cb::pfba(&reference, 1.0).and_then(|w| cb::linear_moma(&perturbed, &w.fluxes)),
            "room" => {
                cb::pfba(&reference, 1.0).and_then(|w| cb::room(&perturbed, &w.fluxes, cb::RoomOptions::default()))
            }
            "loopless" => cb::loopless_fba(&perturbed, 1000.0),
            other => return Err(format!("unknown method {other}")),
        }
        .map_err(text)?;
        let optimal = s.status == Status::Optimal;
        let growth =
            if optimal { Some(perturbed.objective.iter().map(|&(j, c)| c * s.fluxes[j]).sum::<f64>()) } else { None };
        let mut disabled = Vec::new();
        for (r, b) in perturbed.reactions.iter().zip(&reference.reactions) {
            if r.lb == 0.0 && r.ub == 0.0 && (b.lb != 0.0 || b.ub != 0.0) {
                disabled.push(r.id.clone());
            }
        }
        let total_flux = s.fluxes.iter().map(|v| v.abs()).sum();
        Ok(SolveReply {
            status: s.status,
            objective: optimal.then_some(s.objective),
            growth,
            fluxes: if optimal { s.fluxes } else { vec![0.0; perturbed.reactions.len()] },
            disabled,
            milliseconds: started.elapsed().as_secs_f64() * 1000.0,
            total_flux,
        })
    })
    .await
    .map_err(text)?
}

/// Flux range of one reaction at the given share of the optimum.
#[tauri::command]
pub async fn reaction_range(
    state: State<'_, AppState>,
    id: String,
    fraction: f64,
    perturbation: Perturbation,
) -> Reply<(f64, f64)> {
    let base = current(&state)?;
    tauri::async_runtime::spawn_blocking(move || {
        let mut m = apply_medium_and_objective(&base, &perturbation)?;
        apply_knockouts(&mut m, &perturbation)?;
        let j = m.reaction_index(&id).map_err(text)?;
        let r =
            cb::fva(&m, &cb::FvaOptions { fraction_of_optimum: fraction, reactions: Some(vec![j]) }).map_err(text)?;
        Ok((r[0].minimum, r[0].maximum))
    })
    .await
    .map_err(text)?
}

/// Map of the open model: the hand-placed one, one subsystem, or the reactions with most flux.
#[tauri::command]
pub async fn layout(
    state: State<'_, AppState>,
    curated: bool,
    subsystem: Option<String>,
    fluxes: Option<Vec<f64>>,
    count: Option<usize>,
) -> Reply<Value> {
    if curated {
        return serde_json::from_str(CORE_LAYOUT).map_err(text);
    }
    let m = current(&state)?;
    tauri::async_runtime::spawn_blocking(move || {
        let chosen: Vec<usize> = if let Some(sub) = subsystem {
            (0..m.reactions.len()).filter(|&j| m.reactions[j].subsystem.as_deref() == Some(sub.as_str())).collect()
        } else {
            let values = fluxes.unwrap_or_default();
            let mut order: Vec<usize> =
                (0..m.reactions.len()).filter(|&j| !m.reactions[j].id.starts_with("BIOMASS")).collect();
            order.sort_by(|&a, &b| {
                let (x, y) = (values.get(a).copied().unwrap_or(0.0).abs(), values.get(b).copied().unwrap_or(0.0).abs());
                y.partial_cmp(&x).unwrap_or(std::cmp::Ordering::Equal)
            });
            order.truncate(count.unwrap_or(120));
            order
        };
        Ok(computed_layout(&m, &chosen))
    })
    .await
    .map_err(text)?
}

/// A knockout screen.
#[derive(Deserialize)]
pub struct ScreenRequest {
    /// `gene` or `reaction`.
    pub kind: String,
    /// 1 or 2.
    pub order: u8,
    /// Flux whose guaranteed value at optimal growth is reported.
    pub target: Option<String>,
    /// Designs below this objective value are dropped from double screens' starting set.
    pub min_growth: f64,
    #[serde(flatten)]
    pub perturbation: Perturbation,
}

#[derive(Serialize, Clone)]
struct ScreenRow {
    ids: Vec<String>,
    growth: f64,
    guaranteed: Option<f64>,
    possible: Option<f64>,
    disabled: usize,
}

fn screen_sets(m: &Model, req: &ScreenRequest) -> Reply<Vec<(Vec<String>, Vec<usize>)>> {
    // BiGG models name spontaneous reactions with a placeholder that is no real gene.
    let is_real_gene = |id: &str| !matches!(id, "s0001" | "spontaneous" | "SPONTANEOUS" | "S0001");
    let disabled_by_genes = |set: &[usize]| {
        let mut deleted = vec![false; m.genes.len()];
        for &g in set {
            deleted[g] = true;
        }
        m.reactions_disabled_by(&deleted)
    };
    let genes: Vec<usize> = (0..m.genes.len()).filter(|&g| is_real_gene(&m.genes[g].id)).collect();
    let singles_gene = || -> Vec<(Vec<String>, Vec<usize>)> {
        genes.iter().map(|&g| (vec![m.genes[g].id.clone()], disabled_by_genes(&[g]))).collect()
    };
    let singles_reaction = || -> Vec<(Vec<String>, Vec<usize>)> {
        (0..m.reactions.len())
            .filter(|&j| !m.objective.iter().any(|&(o, _)| o == j))
            .map(|j| (vec![m.reactions[j].id.clone()], vec![j]))
            .collect()
    };
    Ok(match (req.kind.as_str(), req.order) {
        ("gene", 1) => singles_gene(),
        ("reaction", 1) => singles_reaction(),
        ("gene", 2) | ("reaction", 2) => {
            let base = if req.kind == "gene" { singles_gene() } else { singles_reaction() };
            // Only sets that disable something can change the answer, and lethal ones cannot recover.
            let viable: Vec<&(Vec<String>, Vec<usize>)> = base
                .iter()
                .filter(|(_, disabled)| !disabled.is_empty())
                .filter(|(_, disabled)| {
                    let mut trial = m.clone();
                    for &j in disabled {
                        trial.reactions[j].lb = 0.0;
                        trial.reactions[j].ub = 0.0;
                    }
                    cb::fba(&trial).is_ok_and(|s| s.status == Status::Optimal && s.objective >= req.min_growth)
                })
                .collect();
            let mut out = Vec::with_capacity(viable.len() * viable.len() / 2);
            for (k, a) in viable.iter().enumerate() {
                for b in &viable[k + 1..] {
                    let mut ids = a.0.clone();
                    ids.extend(b.0.clone());
                    let mut disabled = if req.kind == "gene" {
                        let list: Vec<usize> = ids.iter().filter_map(|g| m.gene_index(g).ok()).collect();
                        disabled_by_genes(&list)
                    } else {
                        let mut d = a.1.clone();
                        d.extend(b.1.iter().copied());
                        d
                    };
                    disabled.sort_unstable();
                    disabled.dedup();
                    out.push((ids, disabled));
                }
            }
            out
        }
        _ => return Err("the kind is gene or reaction and the order is 1 or 2".into()),
    })
}

/// Starts a knockout screen. Each design reports the objective and the target's range at that objective.
#[tauri::command]
pub async fn start_screen(state: State<'_, AppState>, request: ScreenRequest) -> Reply<u64> {
    let base = current(&state)?;
    let jobs = &state.jobs;
    Ok(jobs.spawn(0, move |job| {
        let mut m = apply_medium_and_objective(&base, &request.perturbation)?;
        apply_knockouts(&mut m, &request.perturbation)?;
        let sets = screen_sets(&m, &request)?;
        job.total.store(sets.len(), Ordering::Relaxed);
        let target = match &request.target {
            Some(id) if !id.is_empty() => Some(m.reaction_index(id).map_err(text)?),
            _ => None,
        };
        let mut problem = m.to_lp();
        let objective = problem.obj.clone();
        let sense = problem.sense.unwrap_or(Sense::Maximize);
        let coefs: Vec<(usize, f64)> = objective.iter().enumerate().filter(|(_, c)| **c != 0.0).map(|(j, c)| (j, *c)).collect();
        let retention = problem.add_row(f64::NEG_INFINITY, f64::INFINITY, coefs);
        let backend = lp::default_backend();
        let workers = rayon::current_num_threads().max(1);
        let chunk = (sets.len() / (workers * 2)).max(1);
        let rows: Vec<Result<Vec<ScreenRow>, String>> = sets
            .par_chunks(chunk)
            .map(|slice| {
                let mut session = backend.open(&problem).map_err(text)?;
                let mut out = Vec::with_capacity(slice.len());
                for (ids, disabled) in slice {
                    if job.is_cancelled() {
                        break;
                    }
                    for &j in disabled {
                        session.set_col_bounds(j, 0.0, 0.0).map_err(text)?;
                    }
                    session.set_objective(sense, &objective).map_err(text)?;
                    session.set_row_bounds(retention, f64::NEG_INFINITY, f64::INFINITY).map_err(text)?;
                    let s = session.solve().map_err(text)?;
                    let growth = if s.status == Status::Optimal { s.objective } else { 0.0 };
                    let (mut guaranteed, mut possible) = (None, None);
                    if let (Some(t), Status::Optimal) = (target, s.status) {
                        session.set_row_bounds(retention, growth, f64::INFINITY).map_err(text)?;
                        let mut cost = vec![0.0; objective.len()];
                        cost[t] = 1.0;
                        session.set_objective(Sense::Minimize, &cost).map_err(text)?;
                        let low = session.solve().map_err(text)?;
                        session.set_objective(Sense::Maximize, &cost).map_err(text)?;
                        let high = session.solve().map_err(text)?;
                        if low.status == Status::Optimal {
                            guaranteed = Some(low.objective);
                        }
                        if high.status == Status::Optimal {
                            possible = Some(high.objective);
                        }
                    }
                    for &j in disabled {
                        session.set_col_bounds(j, m.reactions[j].lb, m.reactions[j].ub).map_err(text)?;
                    }
                    out.push(ScreenRow { ids: ids.clone(), growth, guaranteed, possible, disabled: disabled.len() });
                    job.done.fetch_add(1, Ordering::Relaxed);
                }
                Ok(out)
            })
            .collect();
        let mut all = Vec::with_capacity(sets.len());
        for part in rows {
            all.extend(part?);
        }
        let wild = cb::fba(&m).map_err(text)?;
        Ok(json!({"rows": all, "reference": if wild.status == Status::Optimal { Some(wild.objective) } else { None }, "cancelled": job.is_cancelled()}))
    }))
}

/// A sampling request.
#[derive(Deserialize)]
pub struct SamplingRequest {
    pub method: String,
    pub chains: usize,
    pub samples: usize,
    pub thinning: usize,
    pub burn_in: usize,
    pub seed: u64,
    pub fraction: Option<f64>,
    #[serde(flatten)]
    pub perturbation: Perturbation,
}

fn histogram(values: &[f64], bins: usize) -> Value {
    let (lo, hi) = values.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(*v), b.max(*v)));
    let mut counts = vec![0u64; bins];
    if hi > lo {
        for v in values {
            let k = (((v - lo) / (hi - lo)) * bins as f64).floor() as usize;
            counts[k.min(bins - 1)] += 1;
        }
    } else if !values.is_empty() {
        counts[bins / 2] = values.len() as u64;
    }
    json!({"low": lo, "high": hi, "counts": counts})
}

/// Starts flux sampling; the reply holds diagnostics and a histogram per reaction.
#[tauri::command]
pub async fn start_sampling(state: State<'_, AppState>, request: SamplingRequest) -> Reply<u64> {
    let base = current(&state)?;
    Ok(state.jobs.spawn(1, move |job| {
        let mut m = apply_medium_and_objective(&base, &request.perturbation)?;
        apply_knockouts(&mut m, &request.perturbation)?;
        let options = cb::SamplerOptions {
            chains: request.chains,
            samples_per_chain: request.samples,
            thinning: request.thinning,
            burn_in: request.burn_in,
            seed: request.seed,
            fraction_of_optimum: request.fraction,
        };
        let run = match request.method.as_str() {
            "optgp" => cb::optgp(&m, &options),
            _ => cb::chrr(&m, &options),
        }
        .map_err(text)?;
        job.done.store(1, Ordering::Relaxed);
        let histograms: Vec<Value> = (0..m.reactions.len())
            .map(|r| {
                let values: Vec<f64> = run.chains.iter().flat_map(|c| c.iter().map(|s| s[r])).collect();
                histogram(&values, 36)
            })
            .collect();
        Ok(json!({
            "reactions": run.reactions,
            "dimension": run.dimension,
            "diagnostics": run.diagnostics,
            "histograms": histograms,
            "total_samples": request.chains * request.samples,
        }))
    }))
}

/// Progress or result of a job.
#[tauri::command]
pub fn poll_job(state: State<'_, AppState>, id: u64) -> Reply<JobStatus> {
    state.jobs.poll(id).ok_or_else(|| "this analysis is gone".to_string())
}

/// Asks a job to stop.
#[tauri::command]
pub fn cancel_job(state: State<'_, AppState>, id: u64) {
    state.jobs.cancel(id);
}

/// Production envelope.
#[tauri::command]
pub async fn envelope(
    state: State<'_, AppState>,
    control: String,
    target: String,
    points: usize,
    perturbation: Perturbation,
) -> Reply<Vec<cb::EnvelopePoint>> {
    let base = current(&state)?;
    tauri::async_runtime::spawn_blocking(move || {
        let mut m = apply_medium_and_objective(&base, &perturbation)?;
        apply_knockouts(&mut m, &perturbation)?;
        cb::production_envelope(&m, &control, &target, points, None).map_err(text)
    })
    .await
    .map_err(text)?
}

/// Flux variability of all reactions.
#[tauri::command]
pub async fn variability(
    state: State<'_, AppState>,
    fraction: f64,
    perturbation: Perturbation,
) -> Reply<Vec<cb::FvaRange>> {
    let base = current(&state)?;
    tauri::async_runtime::spawn_blocking(move || {
        let mut m = apply_medium_and_objective(&base, &perturbation)?;
        apply_knockouts(&mut m, &perturbation)?;
        cb::fva(&m, &cb::FvaOptions { fraction_of_optimum: fraction, reactions: None }).map_err(text)
    })
    .await
    .map_err(text)?
}

/// Opens a kinetic model: the bundled glycolysis model when no path is given.
#[tauri::command]
pub async fn kinetic_open(state: State<'_, AppState>, path: Option<String>) -> Reply<Value> {
    let model = tauri::async_runtime::spawn_blocking(move || match path {
        Some(p) => KineticModel::read_sbml(p).map_err(text),
        None => KineticModel::from_sbml_str(TEUSINK_XML).map_err(text),
    })
    .await
    .map_err(text)??;
    let parameters: Vec<Value> = model
        .symbols
        .iter()
        .enumerate()
        .filter(|(slot, name)| {
            !name.contains('#')
                && !name.starts_with("__")
                && !model.state_slots.contains(slot)
                && !model.derived.iter().any(|(d, _)| d == slot)
        })
        .filter(|(_, name)| !model.species.iter().any(|s| &s.id == *name))
        .map(|(slot, name)| json!({"name": name, "value": model.initial[slot]}))
        .collect();
    let species: Vec<Value> = model
        .species
        .iter()
        .map(|s| json!({"id": s.id, "initial": model.initial[s.view], "dynamic": model.state_slots.contains(&s.amount) && model.stoich.iter().any(|e| model.state_slots[e.state] == s.amount)}))
        .collect();
    let view = json!({
        "id": model.id,
        "name": model.name,
        "species": species,
        "reactions": model.reactions.iter().map(|r| &r.id).collect::<Vec<_>>(),
        "parameters": parameters,
        "events": model.events.len(),
        "states": model.n_states(),
    });
    *state.kinetic.lock().map_err(text)? = Some(Arc::new(model));
    Ok(view)
}

fn with_overrides(model: &KineticModel, overrides: &[(String, f64)]) -> Reply<KineticModel> {
    let mut m = model.clone();
    for (name, value) in overrides {
        m.set_value(name, *value).map_err(text)?;
    }
    Ok(m)
}

/// Time course of the kinetic model.
#[tauri::command]
pub async fn kinetic_simulate(
    state: State<'_, AppState>,
    end: f64,
    steps: usize,
    method: String,
    overrides: Vec<(String, f64)>,
) -> Reply<Value> {
    let base = current_kinetic(&state)?;
    tauri::async_runtime::spawn_blocking(move || {
        let m = with_overrides(&base, &overrides)?;
        let times: Vec<f64> = (0..=steps).map(|k| end * k as f64 / steps.max(1) as f64).collect();
        let options = SimOptions { method: if method == "bdf" { Method::Bdf } else { Method::Rodas4 }, rtol: 1e-8, atol: 1e-10, ..SimOptions::default() };
        let started = std::time::Instant::now();
        let course = kinetic::simulate(&m, &times, &options).map_err(text)?;
        let series: BTreeMap<&str, Vec<f64>> = m
            .species
            .iter()
            .map(|s| (s.id.as_str(), course.values.iter().map(|v| v[s.view]).collect()))
            .collect();
        let rates: BTreeMap<&str, Vec<f64>> = m
            .reactions
            .iter()
            .map(|r| (r.id.as_str(), course.values.iter().map(|v| v[r.rate]).collect()))
            .collect();
        Ok(json!({
            "time": course.time,
            "series": series,
            "rates": rates,
            "stats": {"steps": course.stats.steps, "rejected": course.stats.rejected, "jacobians": course.stats.jacobians},
            "milliseconds": started.elapsed().as_secs_f64() * 1000.0,
        }))
    })
    .await
    .map_err(text)?
}

fn rows(m: &nalgebra::DMatrix<f64>) -> Vec<Vec<Option<f64>>> {
    (0..m.nrows()).map(|i| (0..m.ncols()).map(|j| m[(i, j)].is_finite().then_some(m[(i, j)])).collect()).collect()
}

/// Steady state and control coefficients of the kinetic model.
#[tauri::command]
pub async fn kinetic_control(state: State<'_, AppState>, overrides: Vec<(String, f64)>) -> Reply<Value> {
    let base = current_kinetic(&state)?;
    tauri::async_runtime::spawn_blocking(move || {
        let m = with_overrides(&base, &overrides)?;
        let steady = kinetic::steady_state(&m, &SteadyOptions { integrate_first: 200.0, ..SteadyOptions::default() })
            .map_err(text)?;
        let c = kinetic::control_analysis(&m, &steady).map_err(text)?;
        let species: BTreeMap<&str, f64> = m.species.iter().map(|s| (s.id.as_str(), steady.values[s.view])).collect();
        let dynamic: Vec<usize> =
            (0..c.states.len()).filter(|&i| (0..c.reactions.len()).any(|j| c.moieties.n[(i, j)] != 0.0)).collect();
        let concentration = nalgebra::DMatrix::from_fn(dynamic.len(), c.reactions.len(), |a, b| {
            c.scaled_concentration_control[(dynamic[a], b)]
        });
        Ok(json!({
            "reactions": c.reactions,
            "states": dynamic.iter().map(|&i| c.states[i].trim_end_matches("#a").to_string()).collect::<Vec<_>>(),
            "fluxes": c.fluxes,
            "species": species,
            "residual": steady.residual,
            "stable": steady.stable,
            "largest_eigenvalue": steady.largest_eigenvalue,
            "flux_control": rows(&c.scaled_flux_control),
            "concentration_control": rows(&concentration),
            "theorems": c.theorems,
            "conserved_moieties": c.states.len() - c.moieties.independent.len(),
        }))
    })
    .await
    .map_err(text)?
}

/// Stochastic ensemble of the kinetic model.
#[tauri::command]
pub async fn kinetic_stochastic(
    state: State<'_, AppState>,
    end: f64,
    steps: usize,
    runs: usize,
    seed: u64,
    epsilon: Option<f64>,
    overrides: Vec<(String, f64)>,
) -> Reply<Value> {
    let base = current_kinetic(&state)?;
    tauri::async_runtime::spawn_blocking(move || {
        let m = with_overrides(&base, &overrides)?;
        let times: Vec<f64> = (0..=steps).map(|k| end * k as f64 / steps.max(1) as f64).collect();
        let algorithm = epsilon.map_or(ssa::Algorithm::Direct, |e| ssa::Algorithm::TauLeap { epsilon: e });
        let e = ssa::ensemble(&m, &times, runs, algorithm, seed).map_err(text)?;
        let pick = |table: &Vec<Vec<f64>>, slot: usize| -> Vec<f64> { table.iter().map(|v| v[slot]).collect() };
        let mean: BTreeMap<&str, Vec<f64>> =
            m.species.iter().map(|s| (s.id.as_str(), pick(&e.mean, s.amount))).collect();
        let sd: BTreeMap<&str, Vec<f64>> = m.species.iter().map(|s| (s.id.as_str(), pick(&e.sd, s.amount))).collect();
        Ok(json!({"time": e.time, "mean": mean, "sd": sd, "runs": e.runs}))
    })
    .await
    .map_err(text)?
}

/// Projects the pFBA fluxes of the open model onto the Biochemistry site's map.
#[tauri::command]
pub async fn atlas(state: State<'_, AppState>, perturbation: Perturbation) -> Reply<Value> {
    let base = current(&state)?;
    tauri::async_runtime::spawn_blocking(move || {
        let mut m = apply_medium_and_objective(&base, &perturbation)?;
        apply_knockouts(&mut m, &perturbation)?;
        let scene: rashevsky::io::overlay::Scene = serde_json::from_str(SCENE).map_err(text)?;
        let rhea: BTreeMap<String, Vec<String>> = serde_json::from_str(RHEA).map_err(text)?;
        let s = cb::pfba(&m, 1.0).map_err(text)?;
        if s.status != Status::Optimal {
            return Err(format!("the model has no flux solution here ({:?})", s.status));
        }
        let objective = m.objective.iter().map(|&(j, c)| c * s.fluxes[j]).sum::<f64>();
        let overlay =
            rashevsky::io::overlay::build_overlay(&m, &s.fluxes, Some(objective), &scene, &rhea).map_err(text)?;
        let raw: Value = serde_json::from_str(SCENE).map_err(text)?;
        Ok(json!({
            "overlay": overlay,
            "scene": {
                "world": raw["world"],
                "domains": raw["domains"],
                "pathways": raw["pathways"],
                "nodes": raw["nodes"],
                "reactions": raw["reactions"],
                "processes": raw["processes"],
                "release": raw["release"],
            },
        }))
    })
    .await
    .map_err(text)?
}

/// Writes fluxes to a CSV or JSON file, by extension.
#[tauri::command]
pub async fn export_fluxes(state: State<'_, AppState>, path: String, fluxes: Vec<f64>) -> Reply<()> {
    let m = current(&state)?;
    let is_json = Path::new(&path).extension().is_some_and(|e| e == "json");
    let body = if is_json {
        serde_json::to_string_pretty(&rashevsky::io::escher::reaction_data(&m, &fluxes)).map_err(text)?
    } else {
        let mut s = String::from("reaction,name,flux\n");
        for (r, v) in m.reactions.iter().zip(&fluxes) {
            s.push_str(&format!("{},\"{}\",{}\n", r.id, r.name.replace('"', "\"\""), v));
        }
        s
    };
    std::fs::write(&path, body).map_err(text)
}

/// Writes an Escher map of the open model with its flux data beside it.
#[tauri::command]
pub async fn export_escher(
    state: State<'_, AppState>,
    path: String,
    fluxes: Vec<f64>,
    subsystem: Option<String>,
) -> Reply<()> {
    let m = current(&state)?;
    tauri::async_runtime::spawn_blocking(move || {
        let chosen: Vec<usize> = (0..m.reactions.len())
            .filter(|&j| !m.reactions[j].id.starts_with("BIOMASS"))
            .filter(|&j| subsystem.as_ref().is_none_or(|s| m.reactions[j].subsystem.as_ref() == Some(s)))
            .collect();
        let layout = rashevsky::io::layout::force_layout(&m, &chosen, 1, 400);
        let doc = rashevsky::io::escher::map_json(&m, &layout, subsystem.as_deref().unwrap_or(&m.id));
        std::fs::write(&path, serde_json::to_string(&doc).map_err(text)?).map_err(text)?;
        let data = PathBuf::from(&path).with_extension("data.json");
        let fluxes = rashevsky::io::escher::reaction_data(&m, &fluxes);
        std::fs::write(data, serde_json::to_string_pretty(&fluxes).map_err(text)?).map_err(text)
    })
    .await
    .map_err(text)?
}

/// Rewrites an SBML file as Level 3 with FBC and groups.
#[tauri::command]
pub async fn export_sbml(input: String, output: String) -> Reply<()> {
    tauri::async_runtime::spawn_blocking(move || {
        let doc = rashevsky::sbml::read::read_file(&input).map_err(text)?;
        std::fs::write(output, rashevsky::sbml::write::write_string(&doc)).map_err(text)
    })
    .await
    .map_err(text)?
}

/// The validation results the release was built with.
#[tauri::command]
pub fn validation() -> Reply<Value> {
    let mut out = serde_json::Map::new();
    for (name, body) in VALIDATION {
        out.insert((*name).to_string(), serde_json::from_str(body).map_err(text)?);
    }
    Ok(Value::Object(out))
}

/// Writes one line to the error stream when the application runs under the test driver.
#[tauri::command]
pub fn qa_log(line: String) {
    if std::env::var_os("RASHEVSKY_QA").is_some() {
        eprintln!("QA {line}");
    }
}
