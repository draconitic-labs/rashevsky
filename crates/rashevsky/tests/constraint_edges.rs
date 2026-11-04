//! Forms of model and bound that the mainstream models do not exercise: exchanges through boundary
//! species, uptake written forwards, negative optima, infinite bounds and minimized objectives.

use rashevsky::Model;
use rashevsky::cb::{self, DfbaOptions, SamplerOptions, Substrate};
use rashevsky::lp::{Sense, Status};
use rashevsky::model::{Metabolite, Reaction};

fn met(id: &str, boundary: bool) -> Metabolite {
    Metabolite { id: id.into(), compartment: "c".into(), boundary, ..Metabolite::default() }
}

fn rxn(id: &str, terms: &[(usize, f64)], lb: f64, ub: f64) -> Reaction {
    Reaction { id: id.into(), lb, ub, stoichiometry: terms.to_vec(), ..Reaction::default() }
}

fn model(mets: Vec<Metabolite>, reactions: Vec<Reaction>, objective: &[(usize, f64)]) -> Model {
    let mut m = Model::new("edge");
    m.metabolites = mets;
    m.reactions = reactions;
    m.objective = objective.to_vec();
    m.sense = Some(Sense::Maximize);
    m.reindex();
    m
}

/// `A_b -> A -> B -> B_b`, exchanges written through boundary species.
fn through_boundary_species() -> Model {
    model(
        vec![met("A_b", true), met("A", false), met("B", false), met("B_b", true)],
        vec![
            rxn("in", &[(0, -1.0), (1, 1.0)], 0.0, 10.0),
            rxn("step", &[(1, -1.0), (2, 1.0)], 0.0, 1000.0),
            rxn("out", &[(2, -1.0), (3, 1.0)], 0.0, 1000.0),
        ],
        &[(2, 1.0)],
    )
}

#[test]
fn exchanges_through_boundary_species_are_not_internal() {
    let m = through_boundary_species();
    assert_eq!(m.exchanges(), vec![0, 2]);
    assert!(!m.is_boundary_reaction(1));
    assert!((cb::fba(&m).unwrap().objective - 10.0).abs() < 1e-9);
    let loopless = cb::loopless_fba(&m, 1000.0).unwrap();
    assert_eq!(loopless.status, Status::Optimal);
    assert!((loopless.objective - 10.0).abs() < 1e-6, "{}", loopless.objective);
}

#[test]
fn a_reaction_between_boundary_species_does_not_stop_loopless() {
    let mut m = through_boundary_species();
    m.reactions.push(rxn("bb", &[(0, -1.0), (3, 1.0)], 0.0, 5.0));
    m.reindex();
    let s = cb::loopless_fba(&m, 1000.0).unwrap();
    assert_eq!(s.status, Status::Optimal);
    assert!((s.objective - 10.0).abs() < 1e-6);
}

#[test]
fn minimal_medium_counts_uptake_through_a_boundary_species() {
    let m = through_boundary_species();
    let medium = cb::minimal_medium(&m, 1.0, false).unwrap();
    assert!((medium.uptake["in"] - 1.0).abs() < 1e-9, "{:?}", medium.uptake);
    assert!((medium.total - 1.0).abs() < 1e-9);
}

#[test]
fn cycle_free_holds_lower_and_upper_bounds() {
    // A loop with a forced minimum keeps that minimum.
    let m = model(
        vec![met("A", false), met("B", false)],
        vec![
            rxn("in", &[(0, 1.0)], 0.0, 10.0),
            rxn("out", &[(0, -1.0)], 0.0, 1000.0),
            rxn("loop1", &[(0, -1.0), (1, 1.0)], 2.0, 100.0),
            rxn("loop2", &[(1, -1.0), (0, 1.0)], 0.0, 100.0),
        ],
        &[(1, 1.0)],
    );
    let plain = cb::fba(&m).unwrap();
    let free = cb::cycle_free(&m, &plain.fluxes).unwrap();
    assert!((free[2] - 2.0).abs() < 1e-9, "{free:?}");
    assert!((free[3] - 2.0).abs() < 1e-9, "{free:?}");
}

#[test]
fn big_m_rows_accept_infinite_bounds() {
    let m = model(
        vec![met("A", false), met("B", false)],
        vec![
            rxn("in", &[(0, 1.0)], 0.0, 10.0),
            rxn("step", &[(0, -1.0), (1, 1.0)], f64::NEG_INFINITY, f64::INFINITY),
            rxn("out", &[(1, -1.0)], 0.0, f64::INFINITY),
        ],
        &[(2, 1.0)],
    );
    assert!((cb::loopless_fba(&m, 1000.0).unwrap().objective - 10.0).abs() < 1e-6);
    let reference = cb::pfba(&m, 1.0).unwrap().fluxes;
    let s = cb::room(&m, &reference, cb::RoomOptions::default()).unwrap();
    assert_eq!(s.status, Status::Optimal);
}

#[test]
fn an_unbounded_side_is_named_rather_than_multiplied() {
    let m = model(
        vec![met("A", false)],
        vec![rxn("in", &[(0, 1.0)], 0.0, f64::INFINITY), rxn("out", &[(0, -1.0)], 0.0, f64::INFINITY)],
        &[(1, 1.0)],
    );
    let error = cb::loopless_fba(&m, 1000.0).unwrap_err().to_string();
    assert!(error.contains("no finite limit"), "{error}");
}

/// Glucose `G` taken up by a reaction with coefficient `sign` and turned into growth with yield 0.1.
fn dfba_model(sign: f64, lb: f64, ub: f64) -> Model {
    model(
        vec![met("G", false)],
        vec![rxn("EX_g", &[(0, sign)], lb, ub), rxn("grow", &[(0, -10.0)], 0.0, 1000.0)],
        &[(1, 1.0)],
    )
}

fn dfba_options(dt: f64) -> DfbaOptions {
    DfbaOptions {
        biomass: 0.1,
        t_end: 20.0,
        dt,
        substrates: vec![Substrate { exchange: "EX_g".into(), initial: 10.0, vmax: 10.0, km: 0.5 }],
        tracked: Vec::new(),
    }
}

#[test]
fn dfba_treats_either_exchange_direction_alike() {
    let usual = cb::dfba(&dfba_model(-1.0, -1000.0, 1000.0), &dfba_options(0.05)).unwrap();
    let forward = cb::dfba(&dfba_model(1.0, -1000.0, 1000.0), &dfba_options(0.05)).unwrap();
    for k in 0..usual.time.len() {
        assert!((usual.biomass[k] - forward.biomass[k]).abs() < 1e-9, "{k}");
        assert!((usual.concentrations[k][0] - forward.concentrations[k][0]).abs() < 1e-9, "{k}");
    }
    // Uptake at vmax gives growth of vmax * yield, not more.
    assert!(usual.growth_rate[0] <= 1.0 + 1e-9, "{}", usual.growth_rate[0]);
}

#[test]
fn dfba_never_exceeds_the_yield_at_any_step_size() {
    // All of the glucose becomes biomass at a yield of 0.1 g per mmol, whatever the step.
    let expected = 0.1 + 10.0 * 0.1;
    for dt in [1.0, 0.25, 0.01] {
        let run = cb::dfba(&dfba_model(-1.0, -1000.0, 1000.0), &dfba_options(dt)).unwrap();
        let end = *run.biomass.last().unwrap();
        assert!((end - expected).abs() < 1e-6 * expected, "dt {dt}: {end} vs {expected}");
        assert!(*run.concentrations.last().unwrap().first().unwrap() < 1e-6);
    }
}

#[test]
fn dfba_keeps_the_models_own_limits() {
    let closed = cb::dfba(&dfba_model(-1.0, 0.0, 1000.0), &dfba_options(0.5)).unwrap();
    assert!(closed.growth_rate.iter().all(|mu| *mu == 0.0), "{:?}", &closed.growth_rate[..3]);
    let limited = cb::dfba(&dfba_model(-1.0, -2.0, 0.0), &dfba_options(0.5)).unwrap();
    assert!((limited.growth_rate[0] - 0.2).abs() < 1e-9, "{}", limited.growth_rate[0]);
}

#[test]
fn dfba_reaches_the_end_time_and_reports_status() {
    let mut options = dfba_options(0.3);
    options.t_end = 1.0;
    let run = cb::dfba(&dfba_model(-1.0, -1000.0, 1000.0), &options).unwrap();
    assert!((run.time.last().unwrap() - 1.0).abs() < 1e-12, "{:?}", run.time);
    assert!(run.status.iter().all(|s| *s == Status::Optimal));
}

#[test]
fn dfba_tracked_species_at_zero_cannot_be_taken_up() {
    let mut options = dfba_options(0.5);
    options.substrates.clear();
    options.tracked = vec![("EX_g".to_string(), 0.0)];
    let run = cb::dfba(&dfba_model(-1.0, -1000.0, 1000.0), &options).unwrap();
    assert!(run.growth_rate.iter().all(|mu| *mu == 0.0), "{:?}", &run.growth_rate[..3]);
}

#[test]
fn dfba_refuses_an_unbounded_problem() {
    // Nothing limits this growth: the exchange is not followed and has no bound.
    let m = model(vec![met("A", false)], vec![rxn("in", &[(0, 1.0)], 0.0, f64::INFINITY), rxn("grow", &[(0, -1.0)], 0.0, f64::INFINITY)], &[(1, 1.0)]);
    let options = DfbaOptions { biomass: 1.0, t_end: 1.0, dt: 0.5, substrates: Vec::new(), tracked: Vec::new() };
    assert!(cb::dfba(&m, &options).unwrap_err().to_string().contains("unbounded"));
}

/// Objective `a + b` with `a` and `b` in `[0, 10]`, not reaching a single reaction's bound.
fn square() -> Model {
    model(
        vec![met("A", false), met("B", false)],
        vec![
            rxn("in_a", &[(0, 1.0)], 0.0, 10.0),
            rxn("in_b", &[(1, 1.0)], 0.0, 10.0),
            rxn("a", &[(0, -1.0)], 0.0, 1000.0),
            rxn("b", &[(1, -1.0)], 0.0, 1000.0),
        ],
        &[(2, 1.0), (3, 1.0)],
    )
}

fn sampler(fraction: Option<f64>) -> SamplerOptions {
    SamplerOptions { chains: 4, samples_per_chain: 500, thinning: 50, burn_in: 500, seed: 3, fraction_of_optimum: fraction }
}

fn objective_values(run: &cb::FluxSamples, m: &Model) -> Vec<f64> {
    let c = m.objective_vector();
    run.chains.iter().flatten().map(|s| s.iter().zip(&c).map(|(v, w)| v * w).sum()).collect()
}

#[test]
fn both_samplers_hold_a_multi_reaction_objective_floor() {
    let m = square();
    for run in [cb::chrr(&m, &sampler(Some(0.5))).unwrap(), cb::optgp(&m, &sampler(Some(0.5))).unwrap()] {
        let values = objective_values(&run, &m);
        let lowest = values.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(lowest >= 10.0 - 1e-6, "{lowest}");
    }
}

#[test]
fn optgp_leaves_the_face_its_warmups_lie_on() {
    // At 90 % of 20, the mean of a + b over the region a + b >= 18 in the square is 18.667.
    let m = square();
    let mut options = sampler(Some(0.9));
    options.chains = 2;
    options.samples_per_chain = 2000;
    let run = cb::optgp(&m, &options).unwrap();
    let values = objective_values(&run, &m);
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    assert!((mean - 18.667).abs() < 0.1, "{mean}");
}

#[test]
fn samplers_keep_a_minimized_objective_under_its_ceiling() {
    let mut m = square();
    m.reactions[2].lb = 2.0;
    m.reactions[3].lb = 2.0;
    m.sense = Some(Sense::Minimize);
    m.reindex();
    // The minimum is 4, so at a fraction of 0.5 the objective may reach 8.
    for run in [cb::chrr(&m, &sampler(Some(0.5))).unwrap(), cb::optgp(&m, &sampler(Some(0.5))).unwrap()] {
        let highest = objective_values(&run, &m).iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(highest <= 8.0 + 1e-6, "{highest}");
    }
}

#[test]
fn sampling_without_a_fraction_leaves_the_objective_free() {
    let m = model(
        vec![met("A", false)],
        vec![rxn("in", &[(0, 1.0)], -5.0, 5.0), rxn("out", &[(0, -1.0)], -5.0, 5.0)],
        &[(1, 1.0)],
    );
    let run = cb::chrr(&m, &sampler(None)).unwrap();
    let lowest = run.chains.iter().flatten().map(|s| s[1]).fold(f64::INFINITY, f64::min);
    assert!(lowest < -4.0, "{lowest}");
}

#[test]
fn null_space_agrees_with_elimination_on_small_coefficients() {
    // A metabolite balanced with coefficients of 1e-7 still ties two reactions together.
    let m = model(
        vec![met("A", false), met("trace", false)],
        vec![
            rxn("in", &[(0, 1.0)], 0.0, 10.0),
            rxn("t1", &[(0, -1.0), (1, 1e-7)], 0.0, 10.0),
            rxn("t2", &[(0, 1.0), (1, -1e-7)], 0.0, 10.0),
            rxn("out", &[(0, -1.0)], 0.0, 10.0),
        ],
        &[(3, 1.0)],
    );
    let run = cb::chrr(&m, &sampler(None)).unwrap();
    for s in run.chains.iter().flatten() {
        assert!((s[1] - s[2]).abs() < 1e-6, "{} {}", s[1], s[2]);
    }
}

#[test]
fn retention_holds_for_negative_and_zero_optima() {
    let m = model(
        vec![met("A", false)],
        vec![rxn("in", &[(0, -1.0)], -10.0, -1.0), rxn("out", &[(0, 1.0)], -1000.0, 1000.0)],
        &[(0, 1.0)],
    );
    // The optimum is -1, and at 0.5 of it the flux may fall to -2.
    let range = &cb::fva(&m, &cb::FvaOptions { fraction_of_optimum: 0.5, reactions: Some(vec![0]) }).unwrap()[0];
    assert!((range.minimum + 2.0).abs() < 1e-9, "{}", range.minimum);
    assert!((range.maximum + 1.0).abs() < 1e-9, "{}", range.maximum);

    let mut m = m;
    m.sense = Some(Sense::Minimize);
    // The optimum is -10, and at 0.5 of it the flux may rise to -5.
    let range = &cb::fva(&m, &cb::FvaOptions { fraction_of_optimum: 0.5, reactions: Some(vec![0]) }).unwrap()[0];
    assert!((range.minimum + 10.0).abs() < 1e-9, "{}", range.minimum);
    assert!((range.maximum + 5.0).abs() < 1e-9, "{}", range.maximum);
    assert!((cb::pfba(&m, 0.5).unwrap().fluxes[0] + 5.0).abs() < 1e-9);

    // A fraction of zero keeps nothing, even at an optimum of zero.
    let free = model(
        vec![met("A", false)],
        vec![rxn("in", &[(0, 1.0)], 0.0, 1000.0), rxn("out", &[(0, -1.0)], 0.0, 1000.0)],
        &[(1, 1.0)],
    );
    let mut free = free;
    free.sense = Some(Sense::Minimize);
    let range = &cb::fva(&free, &cb::FvaOptions { fraction_of_optimum: 0.0, reactions: Some(vec![1]) }).unwrap()[0];
    assert!((range.maximum - 1000.0).abs() < 1e-9, "{}", range.maximum);
}

#[test]
fn deletions_report_unbounded_rather_than_lethal() {
    let m = model(
        vec![met("A", false)],
        vec![rxn("in", &[(0, 1.0)], 0.0, f64::INFINITY), rxn("grow", &[(0, -1.0)], 0.0, f64::INFINITY), rxn("spare", &[(0, -1.0)], 0.0, 1.0)],
        &[(1, 1.0)],
    );
    let rows = cb::single_reaction_deletion(&m, None).unwrap();
    let spare = rows.iter().find(|d| d.ids == ["spare"]).unwrap();
    assert_eq!(spare.status, Status::Unbounded);
    assert_eq!(spare.growth, f64::INFINITY);
}
