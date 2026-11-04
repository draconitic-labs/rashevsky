//! Property-based checks of the optimizers and of SBML round trips.

use proptest::prelude::*;
use rashevsky::Model;
use rashevsky::cb;
use rashevsky::lp::simplex::Simplex;
use rashevsky::lp::{self, Backend, Problem, Sense, Status};

/// A random bounded linear program, always feasible because zero lies inside every row and column.
fn random_lp() -> impl Strategy<Value = Problem> {
    (2usize..9, 1usize..7).prop_flat_map(|(n, m)| {
        (
            prop::collection::vec(-5.0f64..5.0, n),
            prop::collection::vec(prop::collection::vec(-3.0f64..3.0, n), m),
            prop::collection::vec(0.5f64..8.0, m),
            prop::collection::vec(0.5f64..6.0, n),
            any::<bool>(),
        )
            .prop_map(move |(obj, a, rhs, ub, maximize)| {
                let mut p = Problem {
                    sense: Some(if maximize { Sense::Maximize } else { Sense::Minimize }),
                    ..Problem::default()
                };
                for j in 0..n {
                    p.add_col(obj[j], -ub[j], ub[j], false);
                }
                for i in 0..m {
                    let row: Vec<(usize, f64)> = (0..n).map(|j| (j, a[i][j])).collect();
                    p.add_row(-rhs[i], rhs[i], row);
                }
                p
            })
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    /// Two unrelated solvers must find the same optimum, and each solution must be feasible.
    #[test]
    fn simplex_and_highs_agree(p in random_lp()) {
        let mut a = lp::default_backend().open(&p).unwrap();
        let mut b = Simplex::default().open(&p).unwrap();
        let sa = a.solve().unwrap();
        let sb = b.solve().unwrap();
        prop_assert_eq!(sa.status, Status::Optimal);
        prop_assert_eq!(sb.status, Status::Optimal);
        prop_assert!((sa.objective - sb.objective).abs() < 1e-6, "{} {}", sa.objective, sb.objective);
        for (i, row) in p.rows.iter().enumerate() {
            let activity: f64 = row.iter().map(|&(j, c)| c * sb.x[j]).sum();
            prop_assert!(activity >= p.row_lb[i] - 1e-7 && activity <= p.row_ub[i] + 1e-7);
        }
        for j in 0..p.n_cols() {
            prop_assert!(sb.x[j] >= p.col_lb[j] - 1e-7 && sb.x[j] <= p.col_ub[j] + 1e-7);
        }
    }

    /// A warm start after a bound change must reach the same optimum as a cold start.
    #[test]
    fn warm_start_matches_cold_start(p in random_lp(), col in 0usize..8, fraction in 0.0f64..1.0) {
        let col = col % p.n_cols();
        let new_ub = p.col_lb[col] + fraction * (p.col_ub[col] - p.col_lb[col]);
        for backend in [lp::default_backend(), Box::new(Simplex::default()) as Box<dyn Backend>] {
            let mut warm = backend.open(&p).unwrap();
            warm.solve().unwrap();
            warm.set_col_bounds(col, p.col_lb[col], new_ub).unwrap();
            let w = warm.solve().unwrap();
            let mut changed = p.clone();
            changed.col_ub[col] = new_ub;
            let c = backend.open(&changed).unwrap().solve().unwrap();
            prop_assert_eq!(w.status, c.status);
            if w.status == Status::Optimal {
                prop_assert!((w.objective - c.objective).abs() < 1e-6, "{} {}", w.objective, c.objective);
            }
        }
    }

    /// Tightening a bound can never improve a maximum on the core model.
    #[test]
    fn tightening_never_helps(limit in 0.0f64..10.0) {
        let mut m = core();
        let base = cb::fba(&m).unwrap().objective;
        m.set_bounds("EX_glc__D_e", -limit, 1000.0).unwrap();
        let s = cb::fba(&m).unwrap();
        if s.status == Status::Optimal {
            prop_assert!(s.objective <= base + 1e-9);
        }
    }
}

fn core() -> Model {
    Model::read_sbml(concat!(env!("CARGO_MANIFEST_DIR"), "/../../models/e_coli_core/e_coli_core.xml.gz")).unwrap()
}

#[test]
fn sbml_round_trip_keeps_the_model() {
    let doc = rashevsky::sbml::read::read_file(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../models/e_coli_core/e_coli_core.xml.gz"
    ))
    .unwrap();
    let text = rashevsky::sbml::write::write_string(&doc).unwrap();
    let back = Model::from_sbml_str(&text).unwrap();
    let original = Model::from_document(&doc).unwrap();
    assert_eq!(back.reactions.len(), original.reactions.len());
    assert_eq!(back.metabolites.len(), original.metabolites.len());
    assert_eq!(back.genes.len(), original.genes.len());
    for (a, b) in back.reactions.iter().zip(&original.reactions) {
        assert_eq!(a.id, b.id);
        assert_eq!(a.lb, b.lb);
        assert_eq!(a.ub, b.ub);
        assert_eq!(a.stoichiometry, b.stoichiometry);
    }
    let (x, y) = (cb::fba(&back).unwrap().objective, cb::fba(&original).unwrap().objective);
    assert!((x - y).abs() < 1e-9);
    // Every gene deletion gives the same growth after the round trip.
    let d1 = cb::single_gene_deletion(&back, None).unwrap();
    let d2 = cb::single_gene_deletion(&original, None).unwrap();
    for (a, b) in d1.iter().zip(&d2) {
        assert_eq!(a.ids, b.ids);
        assert!((a.growth - b.growth).abs() < 1e-8);
    }
}

#[test]
fn a_tiny_coefficient_still_constrains_the_flux() {
    // A very fast enzyme: e = 5e-11 v with e <= 1e-9, so v <= 20. Dropping the entry would give v = 1000.
    let mut p = Problem { sense: Some(Sense::Maximize), ..Problem::default() };
    let v = p.add_col(1.0, 0.0, 1000.0, false);
    let e = p.add_col(0.0, 0.0, 1e-9, false);
    p.add_row(0.0, 0.0, vec![(e, 1.0), (v, -5e-11)]);
    let mut session = lp::default_backend().open(&p).unwrap();
    let s = session.solve().unwrap();
    assert_eq!(s.status, Status::Optimal);
    assert!((s.objective - 20.0).abs() < 1e-6, "{}", s.objective);
}
