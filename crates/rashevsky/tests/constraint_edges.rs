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
    let m = model(
        vec![met("A", false)],
        vec![rxn("in", &[(0, 1.0)], 0.0, f64::INFINITY), rxn("grow", &[(0, -1.0)], 0.0, f64::INFINITY)],
        &[(1, 1.0)],
    );
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
    SamplerOptions {
        chains: 4,
        samples_per_chain: 500,
        thinning: 50,
        burn_in: 500,
        seed: 3,
        fraction_of_optimum: fraction,
    }
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
        vec![
            rxn("in", &[(0, 1.0)], 0.0, f64::INFINITY),
            rxn("grow", &[(0, -1.0)], 0.0, f64::INFINITY),
            rxn("spare", &[(0, -1.0)], 0.0, 1.0),
        ],
        &[(1, 1.0)],
    );
    let rows = cb::single_reaction_deletion(&m, None).unwrap();
    let spare = rows.iter().find(|d| d.ids == ["spare"]).unwrap();
    assert_eq!(spare.status, Status::Unbounded);
    assert_eq!(spare.growth, f64::INFINITY);
}

fn read(text: &str) -> Model {
    let doc = rashevsky::sbml::read::read_str(text).unwrap();
    Model::from_document(&doc).unwrap()
}

const COBRA_LEVEL2: &str = r#"<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4"><model id="m">
<listOfCompartments><compartment id="c"/></listOfCompartments>
<listOfSpecies><species id="A" compartment="c"/></listOfSpecies>
<listOfReactions>
<reaction id="in" reversible="true"><listOfProducts><speciesReference species="A"/></listOfProducts>
<notes><html xmlns="http://www.w3.org/1999/xhtml"><p>GENE_ASSOCIATION: b1 and b2</p><p>SUBSYSTEM: X</p></html></notes>
<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>F</ci></math><listOfParameters>
<parameter id="LOWER_BOUND" value="-10"/><parameter id="UPPER_BOUND" value="5"/><parameter id="OBJECTIVE_COEFFICIENT" value="0"/>
</listOfParameters></kineticLaw></reaction>
<reaction id="out" reversible="false"><listOfReactants><speciesReference species="A"/></listOfReactants>
<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>F</ci></math><listOfParameters>
<parameter id="LOWER_BOUND" value="0"/><parameter id="UPPER_BOUND" value="1000"/><parameter id="OBJECTIVE_COEFFICIENT" value="1"/>
</listOfParameters></kineticLaw></reaction>
</listOfReactions></model></sbml>"#;

#[test]
fn cobra_level_2_keeps_bounds_objective_and_gene_rules_in_old_places() {
    let m = read(COBRA_LEVEL2);
    assert_eq!((m.reactions[0].lb, m.reactions[0].ub), (-10.0, 5.0));
    assert!((cb::fba(&m).unwrap().objective - 5.0).abs() < 1e-9);
    // Paragraphs on one line stay apart, and the rule is exactly two genes.
    assert_eq!(m.genes.len(), 2, "{:?}", m.genes.iter().map(|g| &g.id).collect::<Vec<_>>());
}

#[test]
fn fbc_version_1_bounds_and_gene_associations_are_read() {
    let text = r#"<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version1" level="3" version="1" fbc:required="false">
<model id="m"><listOfCompartments><compartment id="c" constant="true" size="1"/></listOfCompartments>
<listOfSpecies><species id="A" compartment="c" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/></listOfSpecies>
<listOfReactions>
<reaction id="in" reversible="true" fast="false"><listOfProducts><speciesReference species="A" constant="true" stoichiometry="1"/></listOfProducts></reaction>
<reaction id="out" reversible="false" fast="false"><listOfReactants><speciesReference species="A" constant="true" stoichiometry="1"/></listOfReactants></reaction>
</listOfReactions>
<fbc:listOfFluxBounds>
<fbc:fluxBound fbc:reaction="in" fbc:operation="less" fbc:value="7"/>
<fbc:fluxBound fbc:reaction="in" fbc:operation="greater" fbc:value="-7"/>
</fbc:listOfFluxBounds>
<fbc:listOfGeneAssociations><fbc:geneAssociation fbc:id="g" fbc:reaction="in"><fbc:and><fbc:gene fbc:reference="b1"/><fbc:gene fbc:reference="b2"/></fbc:and></fbc:geneAssociation></fbc:listOfGeneAssociations>
<fbc:listOfObjectives fbc:activeObjective="o"><fbc:objective fbc:id="o" fbc:type="maximize"><fbc:listOfFluxObjectives><fbc:fluxObjective fbc:reaction="out" fbc:coefficient="1"/></fbc:listOfFluxObjectives></fbc:objective></fbc:listOfObjectives>
</model></sbml>"#;
    let m = read(text);
    assert!((cb::fba(&m).unwrap().objective - 7.0).abs() < 1e-9);
    assert_eq!(m.genes.len(), 2);
}

#[test]
fn fbc_version_3_content_the_flux_model_lacks_is_refused() {
    let base = r#"<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version3" level="3" version="1" fbc:required="false">
<model id="m"><listOfCompartments><compartment id="c" constant="true" size="1"/></listOfCompartments>
<listOfSpecies><species id="A" compartment="c" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/></listOfSpecies>
<fbc:listOfObjectives fbc:activeObjective="o"><fbc:objective fbc:id="o" fbc:type="maximize"><fbc:listOfFluxObjectives><fbc:fluxObjective fbc:reaction="r" fbc:coefficient="1" EXTRA/></fbc:listOfFluxObjectives></fbc:objective></fbc:listOfObjectives>
USER
</model></sbml>"#;
    let quadratic = base.replace("EXTRA", r#"fbc:variableType="quadratic""#).replace("USER", "");
    let doc = rashevsky::sbml::read::read_str(&quadratic).unwrap();
    assert!(Model::from_document(&doc).unwrap_err().to_string().contains("quadratic"));
    let user = base.replace("EXTRA", "").replace(
        "USER",
        r#"<fbc:listOfUserDefinedConstraints><fbc:userDefinedConstraint fbc:lowerBound="a" fbc:upperBound="b"/></fbc:listOfUserDefinedConstraints>"#,
    );
    let doc = rashevsky::sbml::read::read_str(&user).unwrap();
    assert!(Model::from_document(&doc).unwrap_err().to_string().contains("user-defined"));
}

#[test]
fn stoichiometry_math_reaches_the_flux_model() {
    let text = r#"<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4"><model id="m">
<listOfCompartments><compartment id="c"/></listOfCompartments>
<listOfSpecies><species id="A" compartment="c"/></listOfSpecies>
<listOfReactions><reaction id="r" reversible="false"><listOfProducts><speciesReference species="A"><stoichiometryMath><math xmlns="http://www.w3.org/1998/Math/MathML"><cn>3</cn></math></stoichiometryMath></speciesReference></listOfProducts></reaction></listOfReactions></model></sbml>"#;
    assert_eq!(read(text).reactions[0].stoichiometry[0].1, 3.0);
}

#[test]
fn gene_rules_read_operators_without_spaces() {
    let mut names = Vec::new();
    let mut intern = |g: &str| {
        names.push(g.to_string());
        names.len() - 1
    };
    rashevsky::gpr::Gpr::parse("a&&b||c", &mut intern).unwrap();
    assert_eq!(names, ["a", "b", "c"]);
}

#[test]
fn a_task_can_name_a_metabolite_on_both_sides_and_a_bad_reaction_fails_alone() {
    use rashevsky::tasks::{Task, check_tasks};
    let m = model(
        vec![Metabolite { name: "A".into(), ..met("A", false) }, Metabolite { name: "B".into(), ..met("B", false) }],
        vec![rxn("make", &[(0, -1.0), (1, 1.0)], 0.0, 100.0)],
        &[],
    );
    let task = |id: &str, equation: Option<&str>, changed: Option<&str>| Task {
        id: id.into(),
        description: String::new(),
        should_fail: false,
        inputs: vec![("A[c]".to_string(), 0.0, 1.0)],
        outputs: vec![("B[c]".to_string(), 0.0, 1.0)],
        equations: equation.map(|e| (e.to_string(), 0.0, 1.0)).into_iter().collect(),
        changed: changed.map(|c| (c.to_string(), 0.0, 1.0)).into_iter().collect(),
    };
    let tasks = [
        task("both_sides", Some("A[c] => A[c] + B[c]"), None),
        task("bad", None, Some("no_such_reaction")),
        task("plain", None, None),
    ];
    let results = check_tasks(&m, &tasks, false).unwrap();
    assert!(results[0].error.is_none(), "{:?}", results[0]);
    assert!(results[1].error.is_some());
    assert!(results[2].error.is_none(), "{:?}", results[2]);
}

#[test]
fn annotations_survive_a_round_trip() {
    let doc =
        rashevsky::sbml::read::read_file(concat!(env!("CARGO_MANIFEST_DIR"), "/../../models/ich360/iCH360.xml.gz"))
            .unwrap();
    let again = rashevsky::sbml::read::read_str(&rashevsky::sbml::write::write_string(&doc).unwrap()).unwrap();
    let count = |d: &rashevsky::sbml::Document| -> usize {
        d.model.reactions.iter().map(|r| r.annotations.len()).sum::<usize>()
            + d.model.species.iter().map(|s| s.annotations.len()).sum::<usize>()
    };
    assert!(count(&doc) > 1000);
    assert_eq!(count(&doc), count(&again));
}

const REVIEW_MODELS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/models");

fn read_file(name: &str) -> Model {
    Model::read_sbml(format!("{REVIEW_MODELS}/{name}.xml")).unwrap()
}

fn followed(exchange: &str, initial: f64, vmax: f64, km: f64, biomass: f64, t_end: f64, dt: f64) -> DfbaOptions {
    DfbaOptions {
        biomass,
        t_end,
        dt,
        substrates: vec![Substrate { exchange: exchange.into(), initial, vmax, km }],
        tracked: Vec::new(),
    }
}

#[test]
fn dfba_stops_forced_uptake_when_the_substrate_is_gone() {
    // 2 A per unit biomass: 10 mmol/L makes 5 more of biomass, in either way of writing the exchange.
    for name in ["dfba_force_neg", "dfba_force_pos"] {
        let m = read_file(name);
        let run = cb::dfba(&m, &followed("EX_A", 10.0, 1000.0, 0.0, 0.1, 5.0, 0.1)).unwrap();
        let end = *run.biomass.last().unwrap();
        assert!((end - 5.1).abs() < 1e-6, "{name}: {end}");
    }
}

#[test]
fn dfba_uses_all_of_the_substrate_in_one_long_step() {
    // Unbounded uptake and a step of one hour: the closed form is 2 x (e^(mu dt) - 1) = S, so mu = 3.93.
    let m = read_file("dfba_free");
    let run = cb::dfba(&m, &followed("EX_A", 10.0, 1e9, 0.0, 0.1, 3.0, 1.0)).unwrap();
    assert!((run.growth_rate[0] - 3.931_825_6).abs() < 1e-5, "{}", run.growth_rate[0]);
    assert!((run.biomass.last().unwrap() - 5.1).abs() < 1e-6, "{:?}", run.biomass);
}

#[test]
fn dfba_reads_growth_from_the_reaction_and_checks_its_input() {
    // Biomass flux bounded at 0.5, with an objective coefficient of 2.
    let m = read_file("dfba_obj2");
    let run = cb::dfba(&m, &followed("EX_A", 1000.0, 1000.0, 0.0, 1.0, 1.0, 0.01)).unwrap();
    assert!((run.biomass.last().unwrap() - 0.5f64.exp()).abs() < 1e-6, "{:?}", run.biomass.last());

    let mut twice = followed("EX_A", 10.0, 1.0, 1.0, 1.0, 1.0, 0.1);
    twice.tracked = vec![("EX_A".to_string(), 10.0)];
    assert!(cb::dfba(&m, &twice).unwrap_err().to_string().contains("twice"));
    assert!(cb::dfba(&m, &followed("EX_A", 10.0, 1.0, 1.0, -1.0, 1.0, 0.1)).is_err());
    assert!(cb::dfba(&m, &followed("EX_A", -10.0, 1.0, 1.0, 1.0, 1.0, 0.1)).is_err());
    assert!(cb::dfba(&m, &followed("EX_A", 10.0, 1.0, -1.0, 1.0, 1.0, 0.1)).is_err());
}

fn one_route() -> Model {
    model(
        vec![met("A", false)],
        vec![rxn("in", &[(0, 1.0)], 0.0, 5.0), rxn("r", &[(0, -1.0)], 0.0, 1000.0)],
        &[(1, 1.0)],
    )
}

fn enzyme_data(catalysis: Vec<cb::enzyme::Catalysis>, enzymes: Vec<cb::enzyme::Enzyme>, pool: f64) -> cb::EnzymeData {
    cb::EnzymeData { enzymes, catalysis, protein_pool: pool }
}

fn enzyme(id: &str, weight: f64) -> cb::enzyme::Enzyme {
    cb::enzyme::Enzyme { id: id.into(), molecular_weight: weight, abundance: None }
}

fn catalysis(enzyme: &str, subunits: Vec<(String, f64)>, kcat: f64) -> cb::enzyme::Catalysis {
    cb::enzyme::Catalysis {
        reaction: "r".into(),
        enzyme: enzyme.into(),
        subunits,
        kcat_forward: kcat,
        kcat_reverse: 0.0,
    }
}

#[test]
fn enzyme_usage_is_what_the_flux_needs_not_what_the_pool_allows() {
    let d = enzyme_data(vec![catalysis("E", Vec::new(), 100.0)], vec![enzyme("E", 1.0)], 100.0);
    let s = cb::enzyme_constrained_fba(&one_route(), &d).unwrap();
    assert!((s.flux.objective - 5.0).abs() < 1e-8);
    assert!((s.usage["E"] - 0.05).abs() < 1e-8, "{}", s.usage["E"]);
    assert!((s.pool_used - 0.05).abs() < 1e-8, "{}", s.pool_used);
}

#[test]
fn enzyme_data_that_cannot_be_meant_is_refused_and_a_repeated_subunit_adds() {
    let m = one_route();
    let single = |copies: f64| vec![("E".to_string(), copies)];
    for bad in [
        enzyme_data(vec![catalysis("", single(-1.0), 100.0)], vec![enzyme("E", 1.0)], 1.0),
        enzyme_data(vec![catalysis("", single(0.0), 100.0)], vec![enzyme("E", 1.0)], 1.0),
        enzyme_data(vec![catalysis("E", Vec::new(), -100.0)], vec![enzyme("E", 1.0)], 1.0),
        enzyme_data(vec![catalysis("E", Vec::new(), 0.0)], vec![enzyme("E", 1.0)], 1.0),
        enzyme_data(vec![catalysis("E", Vec::new(), 100.0)], vec![enzyme("E", 1.0), enzyme("E", 100.0)], 1.0),
    ] {
        assert!(cb::enzyme_constrained_fba(&m, &bad).is_err());
    }
    // Two entries of E are two copies per cycle, at kcat 100 and weight 1: flux = 100 * 0.05 / 2.
    let twice = vec![("E".to_string(), 1.0), ("E".to_string(), 1.0)];
    let d = enzyme_data(vec![catalysis("", twice, 100.0)], vec![enzyme("E", 1.0)], 0.05);
    let s = cb::enzyme_constrained_fba(&m, &d).unwrap();
    assert!((s.flux.objective - 2.5).abs() < 1e-8, "{}", s.flux.objective);
}

fn rewritten(name: &str) -> String {
    let doc = rashevsky::sbml::read::read_file(format!("{REVIEW_MODELS}/{name}.xml")).unwrap();
    rashevsky::sbml::write::write_string(&doc).unwrap()
}

#[test]
fn writing_keeps_kinetic_law_parameters_and_constraints() {
    let again = read(&rewritten("w_cobra"));
    assert!((cb::fba(&again).unwrap().objective - 10.0).abs() < 1e-9);
    assert!(rewritten("w_con").contains("<listOfConstraints>"));
}

#[test]
fn writing_refuses_what_it_would_silently_change() {
    for name in ["w_udc", "w_quad"] {
        let doc = rashevsky::sbml::read::read_file(format!("{REVIEW_MODELS}/{name}.xml")).unwrap();
        assert!(rashevsky::sbml::write::write_string(&doc).is_err(), "{name}");
    }
    let mut doc = rashevsky::sbml::read::read_file(format!("{REVIEW_MODELS}/w_cobra.xml")).unwrap();
    doc.model.reactions[0].gene_rule = Some("(b1 and b2".to_string());
    assert!(rashevsky::sbml::write::write_string(&doc).is_err());
}

#[test]
fn writing_declares_the_version_its_math_needs_and_leaves_out_empty_identifiers() {
    let text = rewritten("w_v2");
    assert!(text.contains("level3/version2/core"), "{}", &text[..200]);
    assert!(!text.contains("fast="));
    let text = rewritten("z_noid");
    assert!(!text.contains("id=\"\""));
}

#[test]
fn namespaces_keep_to_the_element_that_declares_them() {
    // An association under a default FBC namespace is read.
    let m = read_file("z_gene2");
    assert_eq!(m.genes.len(), 1);
    assert_eq!(cb::single_gene_deletion(&m, None).unwrap()[0].growth, 0.0);
    // A prefix bound in an annotation does not reach what follows it.
    let k = rashevsky::kinetic::KineticModel::read_sbml(format!("{REVIEW_MODELS}/z_leak1.xml")).unwrap();
    let run = rashevsky::kinetic::simulate(&k, &[0.0, 1.0], &rashevsky::kinetic::SimOptions::default()).unwrap();
    assert!((run.series(&k, "S").unwrap()[1] - 3.0 * (-2.0f64).exp()).abs() < 1e-6);
    let m = read_file("w_ns5");
    assert!((cb::fba(&m).unwrap().objective - 10.0).abs() < 1e-9);
}

#[test]
fn a_step_that_overshoots_under_the_pivot_tolerance_never_reads_as_infeasible() {
    use rashevsky::lp::simplex::Simplex;
    use rashevsky::lp::{Backend, Problem, Sense, Status};
    // Growth needs an enzyme at 1e-10 per unit, and the enzyme is limited to 1e-9: growth 10.
    let mut problem = Problem { sense: Some(Sense::Maximize), ..Problem::default() };
    let uptake = problem.add_col(0.0, 0.0, 1000.0, false);
    let enzyme = problem.add_col(0.0, 0.0, 1e-9, false);
    let grow = problem.add_col(1.0, 0.0, 1000.0, false);
    problem.add_row(0.0, 0.0, vec![(uptake, 1.0), (grow, -1.0)]);
    problem.add_row(0.0, 0.0, vec![(enzyme, 1.0), (grow, -1e-10)]);
    let highs = rashevsky::lp::default_backend().open(&problem).unwrap().solve().unwrap();
    assert_eq!(highs.status, Status::Optimal);
    assert!((highs.objective - 10.0).abs() < 1e-6, "{}", highs.objective);
    let simplex = Simplex::default().open(&problem).unwrap().solve().unwrap();
    assert!(
        simplex.status == Status::Unknown
            || (simplex.status == Status::Optimal && (simplex.objective - 10.0).abs() < 1e-6),
        "{:?} {}",
        simplex.status,
        simplex.objective
    );
    // No coefficient is that small here, but the basis makes one: 5e-8 over 100. Growth 20.
    let mut scaled = Problem { sense: Some(Sense::Maximize), ..Problem::default() };
    let supply = scaled.add_col(0.0, 0.0, 0.01000001, false);
    let demand = scaled.add_col(0.0, 1.0, 1.0, false);
    let grow = scaled.add_col(1.0, 0.0, 1000.0, false);
    scaled.add_row(0.0, 0.0, vec![(supply, 100.0), (demand, -1.0), (grow, -5e-8)]);
    let highs = rashevsky::lp::default_backend().open(&scaled).unwrap().solve().unwrap();
    assert!((highs.objective - 20.0).abs() < 1e-4, "{}", highs.objective);
    let simplex = Simplex::default().open(&scaled).unwrap().solve().unwrap();
    assert!(
        simplex.status == Status::Unknown
            || (simplex.status == Status::Optimal && (simplex.objective - 20.0).abs() < 1e-4),
        "{:?} {}",
        simplex.status,
        simplex.objective
    );
    // A dead end with no such move left still is a verdict.
    let mut none = Problem { sense: Some(Sense::Maximize), ..Problem::default() };
    let x = none.add_col(1.0, 1.0, 2.0, false);
    none.add_row(0.0, 0.0, vec![(x, 1.0)]);
    assert_eq!(Simplex::default().open(&none).unwrap().solve().unwrap().status, Status::Infeasible);
}

#[test]
fn an_objective_listing_a_reaction_twice_means_the_same_on_both_backends() {
    use rashevsky::lp::simplex::Simplex;
    let mut m = through_boundary_species();
    m.objective = vec![(2, 1.0), (2, 1.0)];
    let highs = cb::pfba(&m, 1.0).unwrap();
    let simplex = cb::pfba_with(&m, 1.0, &Simplex::default()).unwrap();
    assert_eq!(highs.status, Status::Optimal);
    assert!((highs.objective - simplex.objective).abs() < 1e-9, "{} {}", highs.objective, simplex.objective);
    assert!((highs.fluxes[2] - 10.0).abs() < 1e-9);
}

#[test]
fn repeated_coefficients_in_a_row_add_up_after_a_refactor() {
    use rashevsky::lp::simplex::Simplex;
    use rashevsky::lp::{Backend, Problem, Sense, Status};
    let mut problem = Problem { sense: Some(Sense::Maximize), ..Problem::default() };
    problem.add_col(1.0, 0.0, 10.0, false);
    problem.add_row(f64::NEG_INFINITY, 4.0, vec![(0, 1.0), (0, 1.0)]);
    let mut session = Simplex::default().open(&problem).unwrap();
    assert!((session.solve().unwrap().x[0] - 2.0).abs() < 1e-9);
    // A new row makes the session rebuild its inverse with the column already basic.
    session.add_row(f64::NEG_INFINITY, f64::INFINITY, &[(0, 1.0)]).unwrap();
    let again = session.solve().unwrap();
    assert_eq!(again.status, Status::Optimal);
    assert!((again.x[0] - 2.0).abs() < 1e-9, "{}", again.x[0]);
}
