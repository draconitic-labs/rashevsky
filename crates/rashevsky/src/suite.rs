//! Runner for the SBML Test Suite's semantic cases.
//!
//! A case passes when every output value is within the case's tolerance of the reference:
//! `|sim - ref| <= absolute + relative * |ref|`, the comparison the suite prescribes.

use crate::error::{Error, Result};
use crate::kinetic::{KineticModel, Method, SimOptions, simulate};
use serde::{Deserialize, Serialize};
use std::path::Path;

/// Settings file of one case.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Settings {
    /// Start time.
    pub start: f64,
    /// Duration.
    pub duration: f64,
    /// Number of output steps.
    pub steps: usize,
    /// Output symbols.
    pub variables: Vec<String>,
    /// Absolute tolerance.
    pub absolute: f64,
    /// Relative tolerance.
    pub relative: f64,
    /// Variables reported as amounts.
    pub amount: Vec<String>,
    /// Variables reported as concentrations.
    pub concentration: Vec<String>,
}

/// Parses a case's settings file.
pub fn parse_settings(text: &str) -> Settings {
    let mut s = Settings::default();
    let list = |v: &str| -> Vec<String> {
        v.split(',').map(str::trim).filter(|x| !x.is_empty()).map(str::to_string).collect()
    };
    for line in text.lines() {
        let Some((key, value)) = line.split_once(':') else { continue };
        let value = value.trim();
        match key.trim() {
            "start" => s.start = value.parse().unwrap_or(0.0),
            "duration" => s.duration = value.parse().unwrap_or(0.0),
            "steps" => s.steps = value.parse().unwrap_or(0),
            "variables" => s.variables = list(value),
            "absolute" => s.absolute = value.parse().unwrap_or(0.0),
            "relative" => s.relative = value.parse().unwrap_or(0.0),
            "amount" => s.amount = list(value),
            "concentration" => s.concentration = list(value),
            _ => {}
        }
    }
    s
}

/// Reference results as `(headers, rows)`.
pub fn parse_results(text: &str) -> Result<(Vec<String>, Vec<Vec<f64>>)> {
    let mut lines = text.lines();
    let head = lines.next().ok_or_else(|| Error::Invalid("empty results file".into()))?;
    let headers: Vec<String> = head.split(',').map(|s| s.trim().to_string()).collect();
    let mut rows = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let row = line
            .split(',')
            .map(|v| match v.trim() {
                "INF" | "inf" | "Infinity" => Ok(f64::INFINITY),
                "-INF" | "-inf" | "-Infinity" => Ok(f64::NEG_INFINITY),
                "NaN" | "nan" | "NAN" => Ok(f64::NAN),
                t => t.parse::<f64>().map_err(|_| Error::Invalid(format!("bad number `{t}` in results"))),
            })
            .collect::<Result<Vec<f64>>>()?;
        rows.push(row);
    }
    Ok((headers, rows))
}

/// Result of one case.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CaseResult {
    /// Case number.
    pub case: String,
    /// Test type in the case's description.
    pub test_type: String,
    /// Component tags.
    pub component_tags: Vec<String>,
    /// Test tags.
    pub test_tags: Vec<String>,
    /// Whether the case belongs to the engine's scope.
    pub applicable: bool,
    /// Outcome.
    pub passed: bool,
    /// Why it failed, when it did.
    pub detail: Option<String>,
    /// Wall-clock time.
    pub seconds: f64,
}

fn close(sim: f64, reference: f64, absolute: f64, relative: f64) -> bool {
    if reference.is_nan() {
        return sim.is_nan();
    }
    if reference.is_infinite() {
        return sim == reference;
    }
    (sim - reference).abs() <= absolute + relative * reference.abs()
}

fn best_sbml(dir: &Path, case: &str) -> Option<std::path::PathBuf> {
    for version in ["l3v2", "l3v1", "l2v5", "l2v4", "l2v3", "l2v2", "l2v1"] {
        let p = dir.join(format!("{case}-sbml-{version}.xml"));
        if p.exists() {
            return Some(p);
        }
    }
    None
}

fn time_course(path: &Path, settings: &Settings, results: &(Vec<String>, Vec<Vec<f64>>), method: Method) -> Result<Option<String>> {
    let model = KineticModel::read_sbml(path)?;
    let times: Vec<f64> = results.1.iter().map(|r| r[0]).collect();
    let options = SimOptions { method, ..SimOptions::default() };
    let course = simulate(&model, &times, &options)?;
    for (c, name) in results.0.iter().enumerate().skip(1) {
        let series = if settings.amount.iter().any(|s| s == name) {
            match course.series(&model, &format!("{name}#a")) {
                Ok(v) => v,
                Err(_) => course.series(&model, name)?,
            }
        } else if settings.concentration.iter().any(|s| s == name) {
            match course.concentration(&model, name) {
                Ok(v) => v,
                Err(_) => course.series(&model, name)?,
            }
        } else {
            course.series(&model, name)?
        };
        for (k, row) in results.1.iter().enumerate() {
            if !close(series[k], row[c], settings.absolute, settings.relative) {
                return Ok(Some(format!("{name} at t = {}: got {}, expected {}", row[0], series[k], row[c])));
            }
        }
    }
    Ok(None)
}

/// FBC cases: objective and fluxes at the reference tolerance, by the column names of the results.
fn flux_balance(path: &Path, settings: &Settings, results: &(Vec<String>, Vec<Vec<f64>>)) -> Result<Option<String>> {
    let doc = crate::sbml::read::read_file(path)?;
    let model = crate::Model::from_document(&doc)?;
    let solution = crate::cb::fba(&model)?;
    let row = results.1.first().ok_or_else(|| Error::Invalid("no result row".into()))?;
    for (c, name) in results.0.iter().enumerate() {
        let expected = row[c];
        let found = if name == "time" {
            continue;
        } else if doc.model.objectives.iter().any(|o| &o.id == name) {
            if doc.model.active_objective.as_deref() == Some(name.as_str()) || doc.model.objectives.len() == 1 {
                solution.objective
            } else {
                continue;
            }
        } else if let Some(j) = doc.model.reactions.iter().position(|r| &r.id == name) {
            solution.fluxes[j]
        } else {
            continue;
        };
        if !close(found, expected, settings.absolute, settings.relative) {
            return Ok(Some(format!("{name}: got {found}, expected {expected}")));
        }
    }
    Ok(None)
}

/// Runs one case directory.
pub fn run_case(dir: &Path, method: Method) -> Result<CaseResult> {
    let case = dir.file_name().and_then(|s| s.to_str()).unwrap_or_default().to_string();
    let read = |suffix: &str| std::fs::read_to_string(dir.join(format!("{case}-{suffix}")));
    let settings = parse_settings(&read("settings.txt")?);
    let description = read("model.m").unwrap_or_default();
    let field = |key: &str| -> String {
        description
            .lines()
            .find_map(|l| l.trim().strip_prefix(key).map(str::trim).map(str::to_string))
            .unwrap_or_default()
    };
    let tags = |text: String| -> Vec<String> {
        text.split(',').map(str::trim).filter(|s| !s.is_empty()).map(str::to_string).collect()
    };
    let test_type = field("testType:");
    let component_tags = tags(field("componentTags:"));
    let test_tags = tags(field("testTags:"));
    let is_comp = component_tags.iter().chain(&test_tags).any(|t| t.starts_with("comp:"));
    let applicable = !is_comp;
    let mut result = CaseResult {
        case: case.clone(),
        test_type: test_type.clone(),
        component_tags,
        test_tags,
        applicable,
        passed: false,
        detail: None,
        seconds: 0.0,
    };
    if !applicable {
        result.detail = Some("hierarchical model composition (comp package)".into());
        return Ok(result);
    }
    let Some(path) = best_sbml(dir, &case) else {
        result.detail = Some("no SBML file".into());
        return Ok(result);
    };
    let started = std::time::Instant::now();
    let results = parse_results(&read("results.csv")?)?;
    let outcome = if test_type.starts_with("FluxBalance") {
        flux_balance(&path, &settings, &results)
    } else {
        time_course(&path, &settings, &results, method)
    };
    result.seconds = started.elapsed().as_secs_f64();
    match outcome {
        Ok(None) => result.passed = true,
        Ok(Some(detail)) => result.detail = Some(detail),
        Err(e) => result.detail = Some(e.to_string()),
    }
    Ok(result)
}

/// Runs every case under `root` in parallel. `filter` keeps cases by number.
pub fn run_suite(root: &Path, method: Method, filter: Option<&[String]>) -> Result<Vec<CaseResult>> {
    let mut dirs: Vec<std::path::PathBuf> = std::fs::read_dir(root)?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.is_dir() && p.file_name().and_then(|s| s.to_str()).is_some_and(|s| s.chars().all(|c| c.is_ascii_digit())))
        .filter(|p| {
            filter.is_none_or(|list| {
                let name = p.file_name().and_then(|s| s.to_str()).unwrap_or_default();
                list.iter().any(|c| c == name)
            })
        })
        .collect();
    dirs.sort();
    let parts = crate::parallel::map_chunks(&dirs, 1, |slice| run_case(&slice[0], method));
    parts.into_iter().collect()
}
