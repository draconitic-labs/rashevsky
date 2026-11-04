//! Checks on the bundled iCH360 model. Recorded values come from COBRApy on the same file, as compared
//! in full by validation/v13_ich360_cobrapy.py.

use rashevsky::Model;
use rashevsky::cb::{FvaOptions, fba, fva};

/// Growth on 10 mmol/gDW/h of glucose, from COBRApy.
const GROWTH: f64 = 0.8656270306509;

fn core() -> Model {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../models/ich360/iCH360.xml.gz");
    Model::read_sbml(path).unwrap()
}

#[test]
fn reads_the_published_counts() {
    let m = core();
    assert_eq!(m.reactions.len(), 349);
    assert_eq!(m.metabolites.len(), 304);
    assert_eq!(m.genes.len(), 356);
    assert_eq!(m.exchanges().len(), 26);
}

#[test]
fn aerobic_growth_rate() {
    let s = fba(&core()).unwrap();
    assert!((s.objective - GROWTH).abs() < 1e-8, "{}", s.objective);
}

#[test]
fn pure_rust_simplex_agrees_with_highs() {
    let m = core();
    let a = fba(&m).unwrap();
    let b = rashevsky::cb::fba_with(&m, &rashevsky::lp::simplex::Simplex::default()).unwrap();
    assert!((a.objective - b.objective).abs() < 1e-8, "{} {}", a.objective, b.objective);
}

#[test]
fn a_simplex_session_keeps_its_accuracy_over_many_solves() {
    // Flux variability runs one session through every reaction, and the inverse now survives bound and
    // objective changes, so drift would show here.
    let m = core();
    let options = FvaOptions::default();
    let ours = rashevsky::cb::fva_with(&m, &options, &rashevsky::lp::simplex::Simplex::default()).unwrap();
    let highs = fva(&m, &options).unwrap();
    for (j, (a, b)) in ours.iter().zip(&highs).enumerate() {
        let d = (a.minimum - b.minimum).abs().max((a.maximum - b.maximum).abs());
        assert!(d < 1e-6, "{}: {d}", m.reactions[j].id);
    }
}

#[test]
fn simplex_gene_deletions_match_highs() {
    // A deletion that ends in failure counts as lethal, so a false verdict from a kept basis would read
    // as a lethal gene.
    let m = core();
    let ours = rashevsky::cb::single_gene_deletion_with(&m, None, &rashevsky::lp::simplex::Simplex::default()).unwrap();
    let highs = rashevsky::cb::single_gene_deletion(&m, None).unwrap();
    for (a, b) in ours.iter().zip(&highs) {
        assert_eq!(a.status, b.status, "{:?}", a.ids);
        assert!((a.growth - b.growth).abs() < 1e-6, "{:?}: {} {}", a.ids, a.growth, b.growth);
    }
}

#[test]
fn fva_keeps_growth_fixed() {
    let m = core();
    let ranges = fva(&m, &FvaOptions::default()).unwrap();
    let biomass = m.reaction_index("Biomass").unwrap();
    let r = &ranges[biomass];
    assert!((r.minimum - GROWTH).abs() < 1e-6);
    assert!((r.maximum - GROWTH).abs() < 1e-6);
}

#[test]
fn the_analyses_agree_with_recorded_values() {
    use rashevsky::cb;
    let m = core();
    let total: f64 = cb::pfba(&m, 1.0).unwrap().fluxes.iter().map(|v| v.abs()).sum();
    assert!((total - 758.0250933).abs() < 1e-5, "{total}");
    let loopless = cb::loopless_fba(&m, 1000.0).unwrap();
    assert!((loopless.objective - GROWTH).abs() < 1e-6);
    let medium = cb::minimal_medium(&m, 0.1, false).unwrap();
    assert!((medium.total - 6.4475922).abs() < 1e-6);
    let lethal = cb::single_gene_deletion(&m, None).unwrap().iter().filter(|d| d.growth < 1e-6).count();
    assert_eq!(lethal, 114);
}

#[test]
fn deleting_pgi_reroutes_through_the_pentose_phosphate_pathway() {
    use rashevsky::cb;
    let mut m = core();
    let base = cb::pfba(&m, 1.0).unwrap();
    let g6pdh = m.reaction_index("G6PDH2r").unwrap();
    m.set_bounds("PGI", 0.0, 0.0).unwrap();
    let cut = cb::pfba(&m, 1.0).unwrap();
    // Upper glycolysis is closed, so nearly all of the 10 units of glucose must enter the oxidative
    // pentose phosphate pathway, against about a third of it in the unperturbed cell.
    assert!(base.fluxes[g6pdh] < 6.0, "{}", base.fluxes[g6pdh]);
    assert!(cut.fluxes[g6pdh] > 9.0, "{}", cut.fluxes[g6pdh]);
    let biomass = m.reaction_index("Biomass").unwrap();
    assert!((cut.fluxes[biomass] - 0.8552008).abs() < 1e-6, "{}", cut.fluxes[biomass]);
}

#[test]
fn gene_names_lose_the_sbml_prefix() {
    let mut m = core();
    assert!(m.genes.iter().any(|g| g.id == "b4025"));
    assert!(m.genes.iter().all(|g| !g.id.starts_with("G_")));
    // Phosphoglucose isomerase has two genes in this model, so only deleting both stops it.
    let pgi = m.reaction_index("PGI").unwrap();
    assert_eq!(m.reactions[pgi].gpr.to_rule(&|g| m.genes[g].id.clone()), "b4025 or b2843");
    m.set_bounds("PGI", 0.0, 0.0).unwrap();
    assert!((fba(&m).unwrap().objective - 0.8552008).abs() < 1e-6);
}

#[test]
fn loops_are_removed_while_exchanges_stay() {
    use rashevsky::cb;
    let m = core();
    let plain = cb::fba(&m).unwrap();
    let free = cb::cycle_free(&m, &plain.fluxes).unwrap();
    for (j, r) in m.reactions.iter().enumerate() {
        if m.is_boundary_reaction(j) {
            assert!((free[j] - plain.fluxes[j]).abs() < 1e-9, "{}", r.id);
        }
        assert!(free[j].abs() <= plain.fluxes[j].abs() + 1e-9, "{}", r.id);
    }
}

#[test]
fn sampling_converges_for_any_seed() {
    use rashevsky::cb::{SamplerOptions, chrr};
    let m = core();
    for seed in [1, 2, 3] {
        let options = SamplerOptions {
            chains: 4,
            samples_per_chain: 500,
            thinning: 200,
            burn_in: 5000,
            seed,
            fraction_of_optimum: Some(0.9),
        };
        let s = chrr(&m, &options).unwrap();
        assert!(s.diagnostics.max_split_rhat < 1.05, "seed {seed}: {}", s.diagnostics.max_split_rhat);
        assert!(s.diagnostics.min_ess > 200.0, "seed {seed}: {}", s.diagnostics.min_ess);
        // Every sample keeps growth above the floor.
        let biomass = m.reaction_index("Biomass").unwrap();
        for chain in &s.chains {
            for v in chain {
                assert!(v[biomass] >= 0.9 * GROWTH - 1e-6);
            }
        }
    }
}

#[test]
fn dynamic_fba_consumes_substrate_and_grows() {
    use rashevsky::cb::{DfbaOptions, Substrate, dfba};
    let m = core();
    let options = DfbaOptions {
        biomass: 0.01,
        t_end: 12.0,
        dt: 0.02,
        // With a tiny half-saturation constant, uptake stays at its maximum until the glucose is gone.
        substrates: vec![Substrate { exchange: "EX_glc__D_e".into(), initial: 10.0, vmax: 10.0, km: 1e-6 }],
        tracked: vec![],
    };
    let r = dfba(&m, &options).unwrap();
    let end = r.biomass.last().copied().unwrap();
    // Exponential growth at the FBA rate at full uptake: 0.01 exp(0.8656 t).
    assert!((r.biomass[100] - 0.01 * (GROWTH * 2.0).exp()).abs() < 1e-6, "{}", r.biomass[100]);
    assert!(end > 0.5 && end < 1.0, "{end}");
    assert!(r.concentrations.last().unwrap()[0] < 0.1);
    assert!(r.biomass.windows(2).all(|w| w[1] >= w[0] - 1e-12));
}
