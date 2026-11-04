//! The command line from end to end.

use std::process::Command;

const CORE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../models/e_coli_core/e_coli_core.xml.gz");

fn run(args: &[&str]) -> serde_json::Value {
    let output = Command::new(env!("CARGO_BIN_EXE_rashevsky")).args(args).output().unwrap();
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    serde_json::from_slice(&output.stdout).unwrap()
}

#[test]
fn fba_reports_growth() {
    let out = run(&["fba", CORE]);
    assert!((out["objective"].as_f64().unwrap() - 0.8739215).abs() < 1e-6);
    assert_eq!(out["status"], "Optimal");
}

#[test]
fn deletions_and_bounds_take_effect() {
    let out = run(&["fba", CORE, "--ko-reaction", "PGI"]);
    assert!((out["objective"].as_f64().unwrap() - 0.8631).abs() < 1e-3);
    let anaerobic = run(&["fba", CORE, "--bound", "EX_o2_e=0:1000"]);
    assert!((anaerobic["objective"].as_f64().unwrap() - 0.21166).abs() < 1e-4);
    let gene = run(&["fba", CORE, "--ko-gene", "b4025"]);
    assert!((gene["objective"].as_f64().unwrap() - 0.8631).abs() < 1e-3);
}

#[test]
fn simplex_backend_is_selectable() {
    let out = run(&["fba", CORE, "--backend", "simplex"]);
    assert!((out["objective"].as_f64().unwrap() - 0.8739215).abs() < 1e-6);
}

#[test]
fn an_unknown_reaction_is_an_error_with_a_message() {
    let output =
        Command::new(env!("CARGO_BIN_EXE_rashevsky")).args(["fba", CORE, "--ko-reaction", "NOPE"]).output().unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("NOPE"));
}

#[test]
fn info_counts_the_model() {
    let out = run(&["info", CORE]);
    assert_eq!(out["reactions"], 95);
    assert_eq!(out["genes"], 137);
}
