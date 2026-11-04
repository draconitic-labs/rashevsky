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

#[test]
fn the_analyses_agree_with_recorded_values() {
    use rashevsky::cb;
    let m = core();
    // Recorded from the COBRApy comparison in validation/v01_core_cobrapy.py.
    let total: f64 = cb::pfba(&m, 1.0).unwrap().fluxes.iter().map(|v| v.abs()).sum();
    assert!((total - 518.4220855).abs() < 1e-5, "{total}");
    let loopless = cb::loopless_fba(&m, 1000.0).unwrap();
    assert!((loopless.objective - 0.8739215).abs() < 1e-6);
    let medium = cb::minimal_medium(&m, 0.1, false).unwrap();
    assert!((medium.total - 6.7870575).abs() < 1e-6);
    let lethal = cb::single_gene_deletion(&m, None).unwrap().iter().filter(|d| d.growth < 1e-6).count();
    assert_eq!(lethal, 7);
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
    // pentose phosphate pathway, against about half of it in the unperturbed cell.
    assert!(base.fluxes[g6pdh] < 6.0, "{}", base.fluxes[g6pdh]);
    assert!(cut.fluxes[g6pdh] > 9.0, "{}", cut.fluxes[g6pdh]);
    let biomass = m.reaction_index("BIOMASS_Ecoli_core_w_GAM").unwrap();
    assert!((cut.fluxes[biomass] - 0.8631).abs() < 1e-3, "{}", cut.fluxes[biomass]);
}

#[test]
fn loops_are_removed_while_exchanges_stay() {
    use rashevsky::cb;
    let m = core();
    let plain = cb::fba(&m).unwrap();
    let free = cb::cycle_free(&m, &plain.fluxes).unwrap();
    for (j, r) in m.reactions.iter().enumerate() {
        if r.is_boundary() {
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
        let biomass = m.reaction_index("BIOMASS_Ecoli_core_w_GAM").unwrap();
        for chain in &s.chains {
            for v in chain {
                assert!(v[biomass] >= 0.9 * 0.8739215 - 1e-6);
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
    // Exponential growth at the FBA rate at full uptake: 0.01 exp(0.8739 t).
    assert!((r.biomass[100] - 0.01 * (0.8739215f64 * 2.0).exp()).abs() < 1e-6, "{}", r.biomass[100]);
    assert!(end > 0.5 && end < 1.0, "{end}");
    assert!(r.concentrations.last().unwrap()[0] < 0.1);
    assert!(r.biomass.windows(2).all(|w| w[1] >= w[0] - 1e-12));
}
