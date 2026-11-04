//! Python bindings. Results cross as plain dictionaries and lists, which pandas and NumPy take directly.

use engine::cb;
use engine::kinetic::{self, KineticModel, Method, SimOptions, SteadyOptions};
use engine::lp::Sense;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use serde::Serialize;
use serde_json::{Value, json};
use std::collections::BTreeMap;

fn err(e: engine::Error) -> PyErr {
    match e {
        engine::Error::Invalid(_) | engine::Error::Unknown { .. } => PyValueError::new_err(e.to_string()),
        other => PyRuntimeError::new_err(other.to_string()),
    }
}

fn to_python<'py>(py: Python<'py>, value: &Value) -> PyResult<Bound<'py, PyAny>> {
    Ok(match value {
        Value::Null => py.None().into_bound(py),
        Value::Bool(b) => b.into_pyobject(py)?.to_owned().into_any(),
        Value::Number(n) => match n.as_i64() {
            Some(i) => i.into_pyobject(py)?.into_any(),
            None => n.as_f64().unwrap_or(f64::NAN).into_pyobject(py)?.into_any(),
        },
        Value::String(s) => s.into_pyobject(py)?.into_any(),
        Value::Array(items) => {
            let list = PyList::empty(py);
            for item in items {
                list.append(to_python(py, item)?)?;
            }
            list.into_any()
        }
        Value::Object(map) => {
            let dict = PyDict::new(py);
            for (k, v) in map {
                dict.set_item(k, to_python(py, v)?)?;
            }
            dict.into_any()
        }
    })
}

fn result<'py>(py: Python<'py>, value: &impl Serialize) -> PyResult<Bound<'py, PyAny>> {
    let json = serde_json::to_value(value).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    // serde_json writes NaN as null; flux results use `None` for "not computed".
    to_python(py, &json)
}

/// Converts a dict or JSON-compatible Python object through the `json` module.
fn from_python<T: serde::de::DeserializeOwned>(obj: &Bound<'_, PyAny>) -> PyResult<T> {
    let py = obj.py();
    let text: String = py.import("json")?.call_method1("dumps", (obj,))?.extract()?;
    serde_json::from_str(&text).map_err(|e| PyValueError::new_err(e.to_string()))
}

fn fluxes_json(model: &engine::Model, s: &cb::FluxSolution) -> Value {
    let map = |v: &[f64]| -> BTreeMap<&str, f64> {
        v.iter().zip(&model.reactions).map(|(x, r)| (r.id.as_str(), *x)).collect()
    };
    let prices: BTreeMap<&str, f64> =
        s.shadow_prices.iter().zip(&model.metabolites).map(|(x, m)| (m.id.as_str(), *x)).collect();
    json!({
        "status": s.status,
        "objective": s.objective,
        "fluxes": map(&s.fluxes),
        "shadow_prices": prices,
        "reduced_costs": map(&s.reduced_costs),
    })
}

/// A constraint-based model.
#[pyclass(name = "Model", module = "rashevsky")]
#[derive(Clone)]
struct PyModel {
    inner: engine::Model,
}

impl PyModel {
    /// A copy with the given genes and reactions deleted.
    fn with_knockouts(&self, genes: &[String], reactions: &[String]) -> PyResult<engine::Model> {
        let mut m = self.inner.clone();
        for r in reactions {
            m.set_bounds(r, 0.0, 0.0).map_err(err)?;
        }
        if !genes.is_empty() {
            let mut deleted = vec![false; m.genes.len()];
            for g in genes {
                deleted[m.gene_index(g).map_err(err)?] = true;
            }
            for j in m.reactions_disabled_by(&deleted) {
                m.reactions[j].lb = 0.0;
                m.reactions[j].ub = 0.0;
            }
        }
        Ok(m)
    }
}

#[pymethods]
impl PyModel {
    /// Size, compartments and objective.
    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        result(py, &self.inner.summary())
    }

    /// Reaction identifiers.
    #[getter]
    fn reactions(&self) -> Vec<String> {
        self.inner.reactions.iter().map(|r| r.id.clone()).collect()
    }

    /// Metabolite identifiers.
    #[getter]
    fn metabolites(&self) -> Vec<String> {
        self.inner.metabolites.iter().map(|m| m.id.clone()).collect()
    }

    /// Gene identifiers.
    #[getter]
    fn genes(&self) -> Vec<String> {
        self.inner.genes.iter().map(|g| g.id.clone()).collect()
    }

    /// An independent copy.
    fn copy(&self) -> Self {
        self.clone()
    }

    /// Current bounds of a reaction.
    fn bounds(&self, reaction: &str) -> PyResult<(f64, f64)> {
        let j = self.inner.reaction_index(reaction).map_err(err)?;
        Ok((self.inner.reactions[j].lb, self.inner.reactions[j].ub))
    }

    /// Sets the bounds of a reaction.
    fn set_bounds(&mut self, reaction: &str, lb: f64, ub: f64) -> PyResult<()> {
        self.inner.set_bounds(reaction, lb, ub).map_err(err)
    }

    /// Replaces the objective with one reaction.
    #[pyo3(signature = (reaction, maximize=true))]
    fn set_objective(&mut self, reaction: &str, maximize: bool) -> PyResult<()> {
        let sense = if maximize { Sense::Maximize } else { Sense::Minimize };
        self.inner.set_objective(reaction, sense).map_err(err)
    }

    /// The medium as exchange reaction to uptake limit.
    fn medium(&self) -> BTreeMap<String, f64> {
        cb::medium(&self.inner)
    }

    /// Sets the medium; exchanges not listed are closed for uptake.
    fn set_medium(&mut self, uptake: BTreeMap<String, f64>) -> PyResult<()> {
        cb::set_medium(&mut self.inner, &uptake).map_err(err)
    }

    /// Flux balance analysis, optionally after deleting genes or reactions.
    #[pyo3(signature = (knockout_genes=vec![], knockout_reactions=vec![]))]
    fn fba<'py>(
        &self,
        py: Python<'py>,
        knockout_genes: Vec<String>,
        knockout_reactions: Vec<String>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let m = self.with_knockouts(&knockout_genes, &knockout_reactions)?;
        let s = py.detach(|| cb::fba(&m)).map_err(err)?;
        to_python(py, &fluxes_json(&m, &s))
    }

    /// Parsimonious FBA.
    #[pyo3(signature = (fraction=1.0, knockout_genes=vec![], knockout_reactions=vec![]))]
    fn pfba<'py>(
        &self,
        py: Python<'py>,
        fraction: f64,
        knockout_genes: Vec<String>,
        knockout_reactions: Vec<String>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let m = self.with_knockouts(&knockout_genes, &knockout_reactions)?;
        let s = py.detach(|| cb::pfba(&m, fraction)).map_err(err)?;
        to_python(py, &fluxes_json(&m, &s))
    }

    /// Flux variability analysis as a dictionary of `(minimum, maximum)`.
    #[pyo3(signature = (fraction=1.0, reactions=None, knockout_genes=vec![], knockout_reactions=vec![]))]
    fn fva(
        &self,
        py: Python<'_>,
        fraction: f64,
        reactions: Option<Vec<String>>,
        knockout_genes: Vec<String>,
        knockout_reactions: Vec<String>,
    ) -> PyResult<BTreeMap<String, (f64, f64)>> {
        let m = self.with_knockouts(&knockout_genes, &knockout_reactions)?;
        let targets = match reactions {
            None => None,
            Some(list) => {
                Some(list.iter().map(|r| m.reaction_index(r)).collect::<engine::Result<Vec<_>>>().map_err(err)?)
            }
        };
        let options = cb::FvaOptions { fraction_of_optimum: fraction, reactions: targets };
        let ranges = py.detach(|| cb::fva(&m, &options)).map_err(err)?;
        Ok(ranges.into_iter().map(|r| (r.id, (r.minimum, r.maximum))).collect())
    }

    /// Mixed-integer loopless FBA.
    #[pyo3(signature = (big_m=1000.0))]
    fn loopless<'py>(&self, py: Python<'py>, big_m: f64) -> PyResult<Bound<'py, PyAny>> {
        let m = &self.inner;
        let s = py.detach(|| cb::loopless_fba(m, big_m)).map_err(err)?;
        to_python(py, &fluxes_json(m, &s))
    }

    /// MOMA or linear MOMA after a perturbation, relative to the unperturbed pFBA solution.
    #[pyo3(signature = (knockout_genes=vec![], knockout_reactions=vec![], linear=false))]
    fn moma<'py>(
        &self,
        py: Python<'py>,
        knockout_genes: Vec<String>,
        knockout_reactions: Vec<String>,
        linear: bool,
    ) -> PyResult<Bound<'py, PyAny>> {
        let perturbed = self.with_knockouts(&knockout_genes, &knockout_reactions)?;
        let reference = &self.inner;
        let s = py
            .detach(|| {
                let wild = cb::pfba(reference, 1.0)?;
                if linear { cb::linear_moma(&perturbed, &wild.fluxes) } else { cb::moma(&perturbed, &wild.fluxes) }
            })
            .map_err(err)?;
        to_python(py, &fluxes_json(&perturbed, &s))
    }

    /// ROOM after a perturbation, relative to the unperturbed pFBA solution.
    #[pyo3(signature = (knockout_genes=vec![], knockout_reactions=vec![], delta=0.03, epsilon=0.001, relax=false))]
    fn room<'py>(
        &self,
        py: Python<'py>,
        knockout_genes: Vec<String>,
        knockout_reactions: Vec<String>,
        delta: f64,
        epsilon: f64,
        relax: bool,
    ) -> PyResult<Bound<'py, PyAny>> {
        let perturbed = self.with_knockouts(&knockout_genes, &knockout_reactions)?;
        let reference = &self.inner;
        let s = py
            .detach(|| {
                let wild = cb::pfba(reference, 1.0)?;
                cb::room(&perturbed, &wild.fluxes, cb::RoomOptions { delta, epsilon, relax })
            })
            .map_err(err)?;
        to_python(py, &fluxes_json(&perturbed, &s))
    }

    /// Single gene deletions as gene to objective value.
    fn single_gene_deletion(&self, py: Python<'_>) -> PyResult<BTreeMap<String, f64>> {
        let m = &self.inner;
        let d = py.detach(|| cb::single_gene_deletion(m, None)).map_err(err)?;
        Ok(d.into_iter().map(|x| (x.ids[0].clone(), x.growth)).collect())
    }

    /// Single reaction deletions as reaction to objective value.
    fn single_reaction_deletion(&self, py: Python<'_>) -> PyResult<BTreeMap<String, f64>> {
        let m = &self.inner;
        let d = py.detach(|| cb::single_reaction_deletion(m, None)).map_err(err)?;
        Ok(d.into_iter().map(|x| (x.ids[0].clone(), x.growth)).collect())
    }

    /// Double gene deletions as `(gene, gene, objective)`, over the given genes or all of them.
    #[pyo3(signature = (genes=None))]
    fn double_gene_deletion(&self, py: Python<'_>, genes: Option<Vec<String>>) -> PyResult<Vec<(String, String, f64)>> {
        let m = &self.inner;
        let list = match genes {
            None => None,
            Some(ids) => Some(ids.iter().map(|g| m.gene_index(g)).collect::<engine::Result<Vec<_>>>().map_err(err)?),
        };
        let d = py.detach(|| cb::double_gene_deletion(m, list.as_deref())).map_err(err)?;
        Ok(d.into_iter().map(|x| (x.ids[0].clone(), x.ids[1].clone(), x.growth)).collect())
    }

    /// Synthetic lethal gene pairs.
    #[pyo3(signature = (threshold=1e-3))]
    fn synthetic_lethal(&self, py: Python<'_>, threshold: f64) -> PyResult<Vec<(String, String)>> {
        let m = &self.inner;
        let d = py.detach(|| cb::synthetic_lethal_genes(m, threshold)).map_err(err)?;
        Ok(d.into_iter().map(|x| (x.ids[0].clone(), x.ids[1].clone())).collect())
    }

    /// Production envelope of `target` against `control`, as a list of dictionaries.
    #[pyo3(signature = (control, target, points=20))]
    fn production_envelope<'py>(
        &self,
        py: Python<'py>,
        control: &str,
        target: &str,
        points: usize,
    ) -> PyResult<Bound<'py, PyAny>> {
        let m = &self.inner;
        let e = py.detach(|| cb::production_envelope(m, control, target, points, None)).map_err(err)?;
        result(py, &e)
    }

    /// Minimal medium for a growth rate.
    #[pyo3(signature = (growth, minimize_components=false))]
    fn minimal_medium<'py>(
        &self,
        py: Python<'py>,
        growth: f64,
        minimize_components: bool,
    ) -> PyResult<Bound<'py, PyAny>> {
        let m = &self.inner;
        let r = py.detach(|| cb::minimal_medium(m, growth, minimize_components)).map_err(err)?;
        result(py, &r)
    }

    /// Dynamic FBA from a settings dictionary with the fields of `DfbaOptions`.
    fn dfba<'py>(&self, py: Python<'py>, settings: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let options: cb::DfbaOptions = from_python(settings)?;
        let m = &self.inner;
        let r = py.detach(|| cb::dfba(m, &options)).map_err(err)?;
        result(py, &r)
    }

    /// Uniform flux samples. `method` is `chrr` or `optgp`. Diagnostics only, unless `full`.
    #[pyo3(signature = (method="chrr", chains=4, samples=1000, thinning=100, burn_in=2000, seed=1, fraction=None, full=false))]
    #[allow(clippy::too_many_arguments)]
    fn sample<'py>(
        &self,
        py: Python<'py>,
        method: &str,
        chains: usize,
        samples: usize,
        thinning: usize,
        burn_in: usize,
        seed: u64,
        fraction: Option<f64>,
        full: bool,
    ) -> PyResult<Bound<'py, PyAny>> {
        let options = cb::SamplerOptions {
            chains,
            samples_per_chain: samples,
            thinning,
            burn_in,
            seed,
            fraction_of_optimum: fraction,
        };
        let m = &self.inner;
        let run = py
            .detach(|| match method {
                "chrr" => cb::chrr(m, &options),
                "optgp" => cb::optgp(m, &options),
                other => Err(engine::Error::Invalid(format!("unknown sampler `{other}`"))),
            })
            .map_err(err)?;
        if full {
            result(py, &run)
        } else {
            result(py, &json!({"reactions": run.reactions, "dimension": run.dimension, "diagnostics": run.diagnostics}))
        }
    }

    /// Thermodynamics-based flux analysis; `energies` has the fields of `ThermoData`.
    #[pyo3(signature = (energies, c_min=1e-6, c_max=0.02, sigmas=2.0))]
    fn tfa<'py>(
        &self,
        py: Python<'py>,
        energies: &Bound<'py, PyAny>,
        c_min: f64,
        c_max: f64,
        sigmas: f64,
    ) -> PyResult<Bound<'py, PyAny>> {
        let data: cb::ThermoData = from_python(energies)?;
        let options = cb::TfaOptions { c_min, c_max, sigmas, ..cb::TfaOptions::default() };
        let m = &self.inner;
        let s = py.detach(|| cb::tfa(m, &data, &options)).map_err(err)?;
        to_python(
            py,
            &json!({"solution": fluxes_json(m, &s.flux), "concentrations": s.concentrations, "gibbs": s.gibbs}),
        )
    }

    /// Enzyme-constrained FBA; `data` has the fields of `EnzymeData`.
    fn enzyme_fba<'py>(&self, py: Python<'py>, data: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let parsed: cb::EnzymeData = from_python(data)?;
        let m = &self.inner;
        let s = py.detach(|| cb::enzyme_constrained_fba(m, &parsed)).map_err(err)?;
        to_python(py, &json!({"solution": fluxes_json(m, &s.flux), "usage": s.usage, "pool_used": s.pool_used}))
    }

    /// Evaluates a task list in RAVEN format.
    #[pyo3(signature = (path, close_boundaries=true))]
    fn check_tasks<'py>(&self, py: Python<'py>, path: &str, close_boundaries: bool) -> PyResult<Bound<'py, PyAny>> {
        let text = std::fs::read_to_string(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let list = engine::tasks::parse_task_list(&text).map_err(err)?;
        let m = &self.inner;
        let r = py.detach(|| engine::tasks::check_tasks(m, &list, close_boundaries)).map_err(err)?;
        result(py, &r)
    }

    /// Writes the model as SBML Level 3 with FBC.
    fn __repr__(&self) -> String {
        let s = self.inner.summary();
        format!("<Model {}: {} reactions, {} metabolites, {} genes>", s.id, s.reactions, s.metabolites, s.genes)
    }
}

/// Reads an SBML file into a constraint-based model.
#[pyfunction]
fn read_sbml(path: &str) -> PyResult<PyModel> {
    Ok(PyModel { inner: engine::Model::read_sbml(path).map_err(err)? })
}

/// A compiled kinetic model.
#[pyclass(name = "KineticModel", module = "rashevsky")]
#[derive(Clone)]
struct PyKinetic {
    inner: KineticModel,
}

fn method_of(name: &str) -> PyResult<Method> {
    match name {
        "rodas4" => Ok(Method::Rodas4),
        "bdf" => Ok(Method::Bdf),
        other => Err(PyValueError::new_err(format!("unknown method `{other}`"))),
    }
}

#[pymethods]
impl PyKinetic {
    /// Reads and compiles an SBML file.
    #[staticmethod]
    fn read_sbml(path: &str) -> PyResult<Self> {
        Ok(Self { inner: KineticModel::read_sbml(path).map_err(err)? })
    }

    /// Names of all species.
    #[getter]
    fn species(&self) -> Vec<String> {
        self.inner.species.iter().map(|s| s.id.clone()).collect()
    }

    /// Names of all reactions.
    #[getter]
    fn reactions(&self) -> Vec<String> {
        self.inner.reactions.iter().map(|r| r.id.clone()).collect()
    }

    /// Overrides the initial value of a symbol.
    fn set_value(&mut self, name: &str, value: f64) -> PyResult<()> {
        self.inner.set_value(name, value).map_err(err)
    }

    /// Time course over `steps` intervals. `symbols` defaults to all species and reactions.
    #[pyo3(signature = (end, steps=100, method="rodas4", rtol=1e-9, atol=1e-12, symbols=None))]
    #[allow(clippy::too_many_arguments)]
    fn simulate<'py>(
        &self,
        py: Python<'py>,
        end: f64,
        steps: usize,
        method: &str,
        rtol: f64,
        atol: f64,
        symbols: Option<Vec<String>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let options = SimOptions { method: method_of(method)?, rtol, atol, ..SimOptions::default() };
        let times: Vec<f64> = (0..=steps).map(|k| end * k as f64 / steps.max(1) as f64).collect();
        let m = &self.inner;
        let course = py.detach(|| kinetic::simulate(m, &times, &options)).map_err(err)?;
        let names = symbols.unwrap_or_else(|| {
            m.species.iter().map(|s| s.id.clone()).chain(m.reactions.iter().map(|r| r.id.clone())).collect()
        });
        let mut series = BTreeMap::new();
        for name in names {
            let values = course.series(m, &name).map_err(err)?;
            series.insert(name, values);
        }
        result(py, &json!({"time": course.time, "series": series}))
    }

    /// Steady state, with species values and fluxes.
    #[pyo3(signature = (integrate=100.0))]
    fn steady_state<'py>(&self, py: Python<'py>, integrate: f64) -> PyResult<Bound<'py, PyAny>> {
        let m = &self.inner;
        let s = py
            .detach(|| {
                kinetic::steady_state(m, &SteadyOptions { integrate_first: integrate, ..SteadyOptions::default() })
            })
            .map_err(err)?;
        let species: BTreeMap<&str, f64> = m.species.iter().map(|x| (x.id.as_str(), s.values[x.view])).collect();
        let fluxes: BTreeMap<&str, f64> = m.reactions.iter().zip(&s.fluxes).map(|(r, v)| (r.id.as_str(), *v)).collect();
        result(py, &json!({"species": species, "fluxes": fluxes, "residual": s.residual}))
    }

    /// Metabolic control analysis at the steady state; matrices are lists of rows, scaled.
    #[pyo3(signature = (integrate=100.0))]
    fn mca<'py>(&self, py: Python<'py>, integrate: f64) -> PyResult<Bound<'py, PyAny>> {
        let m = &self.inner;
        let c = py
            .detach(|| {
                let state = kinetic::steady_state(
                    m,
                    &SteadyOptions { integrate_first: integrate, ..SteadyOptions::default() },
                )?;
                kinetic::control_analysis(m, &state)
            })
            .map_err(err)?;
        let rows = |a: &rashevsky_matrix::DMatrix<f64>| -> Vec<Vec<f64>> {
            (0..a.nrows()).map(|i| (0..a.ncols()).map(|j| a[(i, j)]).collect()).collect()
        };
        result(
            py,
            &json!({
                "reactions": c.reactions,
                "states": c.states,
                "fluxes": c.fluxes,
                "flux_control": rows(&c.scaled_flux_control),
                "concentration_control": rows(&c.scaled_concentration_control),
                "elasticities": rows(&c.scaled_elasticities),
                "theorems": c.theorems,
            }),
        )
    }

    /// Ensemble of stochastic simulations. `tau_epsilon` selects tau-leaping.
    #[pyo3(signature = (end, steps=100, runs=1000, seed=1, tau_epsilon=None, symbols=None))]
    #[allow(clippy::too_many_arguments)]
    fn ssa<'py>(
        &self,
        py: Python<'py>,
        end: f64,
        steps: usize,
        runs: usize,
        seed: u64,
        tau_epsilon: Option<f64>,
        symbols: Option<Vec<String>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let times: Vec<f64> = (0..=steps).map(|k| end * k as f64 / steps.max(1) as f64).collect();
        let algorithm = match tau_epsilon {
            Some(epsilon) => kinetic::ssa::Algorithm::TauLeap { epsilon },
            None => kinetic::ssa::Algorithm::Direct,
        };
        let m = &self.inner;
        let e = py.detach(|| kinetic::ssa::ensemble(m, &times, runs, algorithm, seed)).map_err(err)?;
        let names = symbols.unwrap_or_else(|| m.species.iter().map(|s| format!("{}#a", s.id)).collect());
        let mut mean = BTreeMap::new();
        let mut sd = BTreeMap::new();
        for name in names {
            let slot = m.slot(&name).map_err(err)?;
            mean.insert(name.clone(), e.mean.iter().map(|v| v[slot]).collect::<Vec<f64>>());
            sd.insert(name, e.sd.iter().map(|v| v[slot]).collect::<Vec<f64>>());
        }
        result(py, &json!({"time": e.time, "mean": mean, "sd": sd, "runs": e.runs}))
    }

    /// Fits parameters. `parameters` lists `{name, start}`, `data` lists `{time, symbol, value, sigma}`.
    fn fit<'py>(
        &self,
        py: Python<'py>,
        parameters: &Bound<'py, PyAny>,
        data: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let parameters: Vec<kinetic::fit::FitParameter> = from_python(parameters)?;
        let data: Vec<kinetic::fit::Observation> = from_python(data)?;
        let m = &self.inner;
        let r = py.detach(|| kinetic::fit::fit(m, &parameters, &data, &SimOptions::default())).map_err(err)?;
        result(py, &r)
    }

    /// Likelihood profile of one fitted parameter, over `factor` either side of its estimate.
    #[pyo3(signature = (parameters, data, parameter, factor=4.0, points=15))]
    fn profile<'py>(
        &self,
        py: Python<'py>,
        parameters: &Bound<'py, PyAny>,
        data: &Bound<'py, PyAny>,
        parameter: &str,
        factor: f64,
        points: usize,
    ) -> PyResult<Bound<'py, PyAny>> {
        let parameters: Vec<kinetic::fit::FitParameter> = from_python(parameters)?;
        let data: Vec<kinetic::fit::Observation> = from_python(data)?;
        let m = &self.inner;
        let options = SimOptions::default();
        let p = py
            .detach(|| {
                let fitted = kinetic::fit::fit(m, &parameters, &data, &options)?;
                kinetic::fit::profile(m, &parameters, &data, &fitted, parameter, factor, points, &options)
            })
            .map_err(err)?;
        result(py, &p)
    }

    fn __repr__(&self) -> String {
        format!(
            "<KineticModel {}: {} species, {} reactions>",
            self.inner.id,
            self.inner.species.len(),
            self.inner.reactions.len()
        )
    }
}

use engine::nalgebra as rashevsky_matrix;

/// Python module.
#[pymodule]
fn rashevsky(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyModel>()?;
    m.add_class::<PyKinetic>()?;
    m.add_function(wrap_pyfunction!(read_sbml, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
