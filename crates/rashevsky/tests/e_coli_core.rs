//! Checks on the bundled `e_coli_core` model.

use rashevsky::Model;
use rashevsky::cb::{FvaOptions, fba, fva};

fn core() -> Model {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../models/e_coli_core/e_coli_core.xml.gz");
    Model::read_sbml(path).unwrap()
}

#[test]
fn reads_bigg_counts() {
    let m = core();
    assert_eq!(m.reactions.len(), 95);
    assert_eq!(m.metabolites.len(), 72);
    assert_eq!(m.genes.len(), 137);
    assert_eq!(m.exchanges().len(), 20);
}

#[test]
fn aerobic_growth_rate() {
    let s = fba(&core()).unwrap();
    assert!((s.objective - 0.8739215069684).abs() < 1e-8, "{}", s.objective);
}

#[test]
fn pure_rust_simplex_agrees_with_highs() {
    let m = core();
    let a = fba(&m).unwrap();
    let b = rashevsky::cb::fba_with(&m, &rashevsky::lp::simplex::Simplex::default()).unwrap();
    assert!((a.objective - b.objective).abs() < 1e-8, "{} {}", a.objective, b.objective);
}

#[test]
fn fva_keeps_growth_fixed() {
    let m = core();
    let ranges = fva(&m, &FvaOptions::default()).unwrap();
    let biomass = m.reaction_index("BIOMASS_Ecoli_core_w_GAM").unwrap();
    let r = &ranges[biomass];
    assert!((r.minimum - 0.8739215069684).abs() < 1e-6);
    assert!((r.maximum - 0.8739215069684).abs() < 1e-6);
}
