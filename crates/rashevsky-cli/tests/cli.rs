//! The command line from end to end.

use std::process::Command;

const CORE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../models/ich360/iCH360.xml.gz");

fn run(args: &[&str]) -> serde_json::Value {
    let output = Command::new(env!("CARGO_BIN_EXE_rashevsky")).args(args).output().unwrap();
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    serde_json::from_slice(&output.stdout).unwrap()
}

#[test]
fn fba_reports_growth() {
    let out = run(&["fba", CORE]);
    assert!((out["objective"].as_f64().unwrap() - 0.8656270).abs() < 1e-6);
    assert_eq!(out["status"], "Optimal");
}

#[test]
fn deletions_and_bounds_take_effect() {
    let out = run(&["fba", CORE, "--ko-reaction", "PGI"]);
    assert!((out["objective"].as_f64().unwrap() - 0.8552008).abs() < 1e-6);
    let anaerobic = run(&["fba", CORE, "--bound", "EX_o2_e=0:1000"]);
    assert!((anaerobic["objective"].as_f64().unwrap() - 0.1565800).abs() < 1e-6);
    // Phosphoglucose isomerase has a second gene, so one deletion changes nothing and both stop it.
    let gene = run(&["fba", CORE, "--ko-gene", "b4025"]);
    assert!((gene["objective"].as_f64().unwrap() - 0.8656270).abs() < 1e-6);
    let both = run(&["fba", CORE, "--ko-gene", "b4025", "--ko-gene", "b2843"]);
    assert!((both["objective"].as_f64().unwrap() - 0.8552008).abs() < 1e-6);
}

#[test]
fn simplex_backend_is_selectable() {
    let out = run(&["fba", CORE, "--backend", "simplex"]);
    assert!((out["objective"].as_f64().unwrap() - 0.8656270).abs() < 1e-6);
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
    assert_eq!(out["reactions"], 349);
    assert_eq!(out["genes"], 356);
}

const BIRTH_DEATH: &str = r#"<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
<model id="birth_death"><listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
<listOfSpecies><species id="S" compartment="c" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/></listOfSpecies>
<listOfParameters><parameter id="k0" value="1.5" constant="true"/><parameter id="k1" value="2" constant="true"/></listOfParameters>
<listOfReactions>
<reaction id="birth" reversible="false"><listOfProducts><speciesReference species="S" stoichiometry="1" constant="true"/></listOfProducts>
<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci>k0</ci><ci>c</ci></apply></math></kineticLaw></reaction>
<reaction id="death" reversible="false"><listOfReactants><speciesReference species="S" stoichiometry="1" constant="true"/></listOfReactants>
<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci>k1</ci><ci>S</ci><ci>c</ci></apply></math></kineticLaw></reaction>
</listOfReactions></model></sbml>"#;

fn scratch(name: &str, text: &str) -> String {
    let dir = std::env::temp_dir().join(format!("rashevsky-cli-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path.to_str().unwrap().to_string()
}

#[test]
fn sensitivity_sobol_and_fit_from_the_command_line() {
    let model = scratch("birth_death.xml", BIRTH_DEATH);
    let local = run(&["sensitivity", &model, "--parameters", "k0,k1"]);
    assert!((local["state"][0][0].as_f64().unwrap() - 0.5).abs() < 1e-8);
    assert!((local["state"][0][1].as_f64().unwrap() + 0.375).abs() < 1e-8);

    let global =
        run(&["sobol", &model, "--range", "k0=1:2", "--range", "k1=1:3:log", "--output", "S", "--samples", "400"]);
    assert_eq!(global["parameters"][1], "k1");
    assert!(global["total"][1].as_f64().unwrap() > global["total"][0].as_f64().unwrap());

    // S(t) from zero is (k0 / k1)(1 - exp(-k1 t)); these rows hold k1 = 2 and k0 = 1.5 exactly.
    let mut table = String::from("time,symbol,value,sigma\n");
    for t in [0.25f64, 0.5, 1.0, 2.0, 4.0] {
        table += &format!("{t},S,{},0.01\n", 0.75 * (1.0 - (-2.0 * t).exp()));
    }
    let data = scratch("measured.csv", &table);
    let fit = run(&["fit", &model, "--start", "k0=4", "--start", "k1=0.5", "--data", &data, "--profile", "k1"]);
    assert!((fit["fit"]["estimates"][0].as_f64().unwrap() - 1.5).abs() < 1e-3, "{fit}");
    assert!((fit["fit"]["estimates"][1].as_f64().unwrap() - 2.0).abs() < 1e-3, "{fit}");
    let (low, high) = (fit["profiles"][0]["lower"].as_f64().unwrap(), fit["profiles"][0]["upper"].as_f64().unwrap());
    assert!(low < 2.0 && high > 2.0, "{fit}");
}

#[test]
fn malformed_ranges_and_tables_are_errors_with_a_message() {
    let model = scratch("birth_death_bad.xml", BIRTH_DEATH);
    let bad_range = Command::new(env!("CARGO_BIN_EXE_rashevsky"))
        .args(["sobol", &model, "--range", "k0=1", "--output", "S"])
        .output()
        .unwrap();
    assert!(!bad_range.status.success());
    assert!(String::from_utf8_lossy(&bad_range.stderr).contains("name=low:high"));
    let table = scratch("bad.csv", "t,s,v\n1,S,2\n");
    let bad_table = Command::new(env!("CARGO_BIN_EXE_rashevsky"))
        .args(["fit", &model, "--start", "k0=1", "--data", &table])
        .output()
        .unwrap();
    assert!(!bad_table.status.success());
    assert!(String::from_utf8_lossy(&bad_table.stderr).contains("time,symbol,value,sigma"));
}
