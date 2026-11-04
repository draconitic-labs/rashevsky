//! Kinetic behaviour outside the mainstream of published models: stochastic simulation of what a
//! propensity cannot describe, steady states of time-dependent and unstable systems, overrides, jumps in
//! rates at or between output times, and delay across events.

use rashevsky::kinetic::fit::{FitParameter, Observation, fit, profile};
use rashevsky::kinetic::sensitivity::local_steady;
use rashevsky::kinetic::ssa::{Algorithm, ensemble};
use rashevsky::kinetic::{KineticModel, Method, SimOptions, SteadyOptions, control_analysis, simulate, steady_state};

const MODELS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/models");

fn model(name: &str) -> KineticModel {
    KineticModel::read_sbml(format!("{MODELS}/{name}.xml")).unwrap()
}

fn end_value(m: &KineticModel, symbol: &str, end: f64, steps: usize, method: Method) -> f64 {
    let times: Vec<f64> = (0..=steps).map(|k| end * k as f64 / steps as f64).collect();
    let options = SimOptions { method, ..SimOptions::default() };
    *simulate(m, &times, &options).unwrap().series(m, symbol).unwrap().last().unwrap()
}

fn mean_at_end(
    m: &KineticModel,
    symbol: &str,
    end: f64,
    steps: usize,
    runs: usize,
    algorithm: Algorithm,
) -> (f64, f64) {
    let times: Vec<f64> = (0..=steps).map(|k| end * k as f64 / steps as f64).collect();
    let e = ensemble(m, &times, runs, algorithm, 11).unwrap();
    let slot = m.slot(symbol).unwrap();
    (e.mean[steps][slot], e.sd[steps][slot])
}

const TAU: Algorithm = Algorithm::TauLeap { epsilon: 0.03 };

#[test]
fn a_leap_that_meets_a_stop_still_fires() {
    // The explicit method's own error at epsilon 0.03 is about 3 %, and it does not grow with stops.
    let exact = 1e6 * (-2.0f64).exp();
    let m = model("tau_decay");
    for steps in [1, 10, 100, 400] {
        let (mean, _) = mean_at_end(&m, "X", 2.0, steps, 40, TAU);
        assert!((mean - exact).abs() < 0.04 * exact, "{steps} stops: {mean} vs {exact}");
    }
    let (mean, _) = mean_at_end(&model("tau_event"), "X", 2.0, 10, 40, TAU);
    assert!((mean - exact / 10.0).abs() < 0.04 * exact / 10.0, "with an event: {mean}");
}

#[test]
fn a_rejected_leap_does_not_advance_the_time() {
    let (mean, sd) = mean_at_end(&model("tau_neg"), "X", 1.0, 2, 400, TAU);
    let exact = 1e5 * (-1.0f64).exp();
    assert!((mean - exact).abs() < 0.01 * exact, "{mean}");
    assert!((sd - (exact * (1.0 - (-1.0f64).exp())).sqrt()).abs() < 30.0, "{sd}");
}

#[test]
fn a_negative_parameter_with_a_rate_rule_does_not_stall_a_leap() {
    mean_at_end(&model("tau_negparam"), "X", 100.0, 1, 5, TAU);
}

#[test]
fn stochastic_simulation_refuses_what_a_propensity_cannot_describe() {
    for (name, expected) in
        [("ssa_rev", "propensity"), ("ssa_delay", "delay"), ("ssa_alg", "algebraic"), ("m_casc", "settle")]
    {
        let m = model(name);
        let error = ensemble(&m, &[0.0, 1.0], 2, Algorithm::Direct, 1).unwrap_err().to_string();
        assert!(error.contains(expected), "{name}: {error}");
    }
}

#[test]
fn a_reaction_cannot_take_more_than_is_there() {
    let m = model("ssa_dimer2");
    let (x, _) = mean_at_end(&m, "X", 100.0, 1, 20, Algorithm::Direct);
    let (d, _) = mean_at_end(&m, "D", 100.0, 1, 20, Algorithm::Direct);
    assert_eq!((x, d), (1.0, 2.0));
    assert_eq!(mean_at_end(&model("ssa_zero"), "X", 100.0, 1, 50, Algorithm::Direct).0, 0.0);
    assert_eq!(mean_at_end(&model("ssa_dimer"), "X", 0.5, 1, 5, Algorithm::Direct).0, 1.0);
}

#[test]
fn stochastic_rates_follow_the_time_and_varying_stoichiometry() {
    let (mean, sd) = mean_at_end(&model("ssa_time"), "X", 100.0, 1, 200, Algorithm::Direct);
    assert!((mean - 500.0).abs() < 5.0 * sd / 200f64.sqrt() + 1.0, "{mean}");
    let (mean, sd) = mean_at_end(&model("ssa_stoich"), "X", 10.0, 1, 200, Algorithm::Direct);
    assert!((mean - 600.0).abs() < 5.0 * sd / 200f64.sqrt() + 1.0, "{mean}");
}

#[test]
fn stochastic_events_respect_priority() {
    assert_eq!(mean_at_end(&model("ssa_prio"), "p", 1.0, 1, 5, Algorithm::Direct).0, 1.0);
}

#[test]
fn the_first_output_time_does_not_shift_the_dynamics() {
    let m = KineticModel::from_sbml_str(DECAY).unwrap();
    let run = simulate(&m, &[5.0, 6.0], &SimOptions::default()).unwrap();
    let a = run.series(&m, "y").unwrap();
    assert!((a[0] - (-5.0f64).exp()).abs() < 1e-7, "{}", a[0]);
    assert!((a[1] - (-6.0f64).exp()).abs() < 1e-7, "{}", a[1]);
    let e = ensemble(&m, &[5.0], 400, Algorithm::Direct, 3).unwrap();
    assert!(e.mean[0][m.slot("y").unwrap()].is_finite());
    assert!(simulate(&m, &[-1.0, 1.0], &SimOptions::default()).is_err());
}

const DECAY: &str = r#"<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2"><model id="m">
<listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
<listOfSpecies><species id="y" compartment="c" initialAmount="1" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/></listOfSpecies>
<listOfReactions><reaction id="r" reversible="false"><listOfReactants><speciesReference species="y" stoichiometry="1" constant="true"/></listOfReactants>
<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>y</ci></math></kineticLaw></reaction></listOfReactions></model></sbml>"#;

#[test]
fn steady_states_read_the_settled_time_and_the_present_delay() {
    let options = SteadyOptions::default();
    let value = |name: &str| {
        let m = model(name);
        let s = steady_state(&m, &options).unwrap();
        s.values[m.slot("S").unwrap()]
    };
    assert!(value("ss_pulse").abs() < 1e-9);
    assert!((value("ss_decay") - 1.0).abs() < 1e-9);
    assert!((value("ss_delay") - 1.0).abs() < 1e-9);
}

#[test]
fn newton_does_not_settle_on_an_unstable_root_when_an_attractor_is_near() {
    let m = model("ss_bist");
    let s = steady_state(&m, &SteadyOptions::default()).unwrap();
    assert!(s.values[m.slot("S").unwrap()].abs() < 1e-9, "{}", s.values[m.slot("S").unwrap()]);
    assert!(s.stable);
}

#[test]
fn a_limit_cycle_is_reported_as_unstable() {
    let s = steady_state(&model("ss_bruss"), &SteadyOptions::default()).unwrap();
    assert!(!s.stable);
    assert!((s.largest_eigenvalue - 0.5).abs() < 1e-6, "{}", s.largest_eigenvalue);
}

#[test]
fn control_analysis_refuses_varying_stoichiometry() {
    let m = model("ss_vs");
    assert!(steady_state(&m, &SteadyOptions::default()).unwrap_err().to_string().contains("stoichiometry"));
}

fn sensitivities(name: &str, parameters: &[&str]) -> rashevsky::kinetic::sensitivity::LocalSensitivity {
    let m = model(name);
    let state = steady_state(&m, &SteadyOptions::default()).unwrap();
    let control = control_analysis(&m, &state).unwrap();
    let names: Vec<String> = parameters.iter().map(|p| (*p).to_string()).collect();
    local_steady(&m, &state, &control, &names).unwrap()
}

#[test]
fn local_sensitivities_reach_local_parameters_and_constant_species() {
    // S = kin X0 / (kloc 2 kg), so at the model's values dS/dkloc = -2 and dS/dX0 = 0.5.
    let local = sensitivities("l_sens", &["out.kloc", "X0", "kin"]);
    assert!((local.state[(0, 0)] + 2.0).abs() < 1e-7, "{}", local.state[(0, 0)]);
    assert!((local.state[(0, 1)] - 0.5).abs() < 1e-7, "{}", local.state[(0, 1)]);
    assert!((local.state[(0, 2)] - 1.0).abs() < 1e-7, "{}", local.state[(0, 2)]);
}

#[test]
fn local_sensitivities_include_a_parameter_in_a_stoichiometry() {
    let local = sensitivities("l_stoich", &["n"]);
    assert!((local.state[(0, 0)] - 1.0).abs() < 1e-7, "{}", local.state[(0, 0)]);
}

#[test]
fn local_sensitivities_refuse_names_they_cannot_vary() {
    let m = model("l_sens");
    let state = steady_state(&m, &SteadyOptions::default()).unwrap();
    let control = control_analysis(&m, &state).unwrap();
    for name in ["nonexistent", "S"] {
        assert!(local_steady(&m, &state, &control, &[name.to_string()]).is_err(), "{name}");
    }
    let rule_set = local_steady(&m, &state, &control, &["kd".to_string()]);
    assert!(rule_set.is_err());
}

#[test]
fn overrides_refuse_rules_and_ignore_their_order() {
    assert!(model("s_rule").set_value("k", 5.0).is_err());
    let mut a = model("s_comp");
    a.set_value("B", 3.0).unwrap();
    a.set_value("c", 4.0).unwrap();
    let mut b = model("s_comp");
    b.set_value("c", 4.0).unwrap();
    b.set_value("B", 3.0).unwrap();
    assert_eq!(a.initial[a.slot("B").unwrap()], 3.0);
    assert_eq!(b.initial[b.slot("B").unwrap()], 3.0);
    assert_eq!(a.initial[a.slot("B#a").unwrap()], 12.0);
    assert_eq!(b.initial[b.slot("B#a").unwrap()], 12.0);
}

#[test]
fn overrides_reach_algebraic_rules_initial_assignments_and_delay() {
    let mut m = model("s_alg");
    m.set_value("k", 3.0).unwrap();
    assert!((m.initial[m.slot("z").unwrap()] - 30.0).abs() < 1e-9);
    let mut m = model("s_ia");
    m.set_value("B#a", 40.0).unwrap();
    assert!((m.initial[m.slot("w").unwrap()] - 20.0).abs() < 1e-9);
    let mut m = model("s_delay");
    m.set_value("k", 1.0).unwrap();
    assert!(end_value(&m, "S", 1.0, 1, Method::Rodas4).is_finite());
}

#[test]
fn an_initial_assignment_reads_what_an_algebraic_rule_solved() {
    let m = model("e_algia");
    assert!((m.initial[m.slot("x").unwrap()] - 2.0).abs() < 1e-9);
}

#[test]
fn delay_history_holds_the_value_before_an_event() {
    // x' = 1, x := 0 at t = 4.8, y = delay(x, 1): y(t) = t - 1 until t = 5.8.
    let m = model("d_event2");
    for method in [Method::Rodas4, Method::Bdf] {
        for (end, expected) in [(5.7, 4.7), (5.79, 4.79)] {
            let y = end_value(&m, "y", end, 579, method);
            assert!((y - expected).abs() < 1e-4, "{method:?} t = {end}: {y}");
        }
    }
}

#[test]
fn delay_differential_equations_keep_their_order() {
    // y' = -y(t - 1) with y = 1 before the start has y(5) = -11/120 + ... by the method of steps.
    let m = model("d_dde");
    let options = SimOptions { rtol: 1e-9, atol: 1e-12, ..SimOptions::default() };
    let run = simulate(&m, &[0.0, 5.0], &options).unwrap();
    assert!(run.stats.steps < 5_000, "{} steps", run.stats.steps);
}

#[test]
fn a_rate_that_jumps_on_or_between_output_times_is_integrated() {
    for method in [Method::Rodas4, Method::Bdf] {
        let label = format!("{method:?}");
        assert!((end_value(&model("sw1"), "S", 20.0, 20, method) - 10.0).abs() < 1e-6, "{label}");
        assert!(end_value(&model("sw1"), "S", 10.0, 10, method).abs() < 1e-6, "{label}");
        assert!((end_value(&model("sw3"), "S", 2000.0, 1, method) - 999.97).abs() < 1e-5, "{label}");
        assert!((end_value(&model("sw4"), "S", 20.0, 20, method) - 10.0).abs() < 1e-6, "{label}");
        assert!(end_value(&model("sw5"), "S", 100.0, 100, method).abs() < 1e-6, "{label}");
    }
}

#[test]
fn a_narrow_window_cannot_fall_between_steps() {
    for method in [Method::Rodas4, Method::Bdf] {
        assert_eq!(end_value(&model("e_window"), "f", 10.0, 1, method), 1.0, "{method:?}");
        assert!((end_value(&model("m_pulse"), "S", 20.0, 1, method) - 0.1).abs() < 1e-6, "{method:?}");
    }
}

#[test]
fn a_compartment_changed_through_a_rule_behaves_as_one_changed_directly() {
    let direct = end_value(&model("e_vdirect"), "S", 2.0, 2, Method::Rodas4);
    let through_rule = end_value(&model("e_vrule"), "S", 2.0, 2, Method::Rodas4);
    assert!((direct - through_rule).abs() < 1e-9, "{direct} {through_rule}");
    assert!((direct - 1.0).abs() < 1e-9, "{direct}");
}

#[test]
fn an_odd_root_of_a_negative_value_has_a_derivative() {
    for method in [Method::Rodas4, Method::Bdf] {
        assert!((end_value(&model("m_root"), "x", 1.0, 1, method) + 10.08115).abs() < 1e-3, "{method:?}");
    }
}

#[test]
fn a_negative_event_delay_is_refused() {
    let m = model("m_negd");
    assert!(simulate(&m, &[0.0, 2.0], &SimOptions::default()).unwrap_err().to_string().contains("delay"));
}

#[test]
fn a_fast_reaction_equilibrates_whatever_the_size_of_its_law() {
    let m = model("fast_slowlaw");
    let run = simulate(&m, &[0.0], &SimOptions::default()).unwrap();
    assert!((run.series(&m, "A").unwrap()[0] - 5.0).abs() < 1e-3, "{:?}", run.series(&m, "A").unwrap());
}

fn observations() -> Vec<Observation> {
    [0.5, 1.0, 2.0, 3.0].iter().map(|&t| Observation { time: t, symbol: "S".into(), value: 0.0, sigma: 0.02 }).collect()
}

#[test]
fn profile_limits_come_from_a_walk_starting_at_the_estimate() {
    let m = model("f_base");
    let exact = |t: f64| 2.0 * (-0.7f64 * t).exp();
    let data: Vec<Observation> = [0.5, 1.0, 2.0, 3.0]
        .iter()
        .map(|&t| Observation { time: t, symbol: "A".into(), value: exact(t), sigma: 0.02 })
        .collect();
    let parameters = [FitParameter { name: "k".into(), start: 1.0 }];
    let options = SimOptions::default();
    let fitted = fit(&m, &parameters, &data, &options).unwrap();
    // With only four grid values, none near the estimate, the limits must still bracket it or be open.
    let p = profile(&m, &parameters, &data, &fitted, "k", 10.0, 4, &options).unwrap();
    assert!(p.values.iter().any(|v| (v - fitted.estimates[0]).abs() < 1e-9));
    if let Some(low) = p.lower {
        assert!(low < fitted.estimates[0]);
    }
    if let Some(high) = p.upper {
        assert!(high > fitted.estimates[0]);
    }
}

#[test]
fn a_failed_simulation_is_not_a_profile_limit() {
    let m = model("f_fail");
    let data = vec![Observation { time: 1.0, symbol: "A".into(), value: 1.0, sigma: 10.0 }];
    let parameters = [FitParameter { name: "k".into(), start: 1.0 }];
    let options = SimOptions::default();
    let fitted = fit(&m, &parameters, &data, &options).unwrap();
    let p = profile(&m, &parameters, &data, &fitted, "k", 4.0, 15, &options).unwrap();
    assert_eq!(p.upper, None, "{:?} {:?}", p.values, p.chi2);
}

#[test]
fn unusable_measurements_and_rule_parameters_are_refused() {
    let m = model("f_base");
    let parameters = [FitParameter { name: "k".into(), start: 1.0 }];
    let nan = vec![Observation { time: 1.0, symbol: "A".into(), value: f64::NAN, sigma: 1.0 }];
    assert!(fit(&m, &parameters, &nan, &SimOptions::default()).is_err());
    let rule = model("f_rule");
    let data = observations();
    assert!(fit(&rule, &parameters, &data, &SimOptions::default()).unwrap_err().to_string().contains("rule"));
}

#[test]
fn sbml_with_a_prefixed_math_namespace_is_read() {
    let m = model("io_prefix");
    assert!((end_value(&m, "A", 1.0, 1, Method::Rodas4) - 2.0 * (-1.0f64).exp()).abs() < 1e-6);
}

#[test]
fn a_required_package_in_use_is_refused() {
    let text = r#"<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:comp="http://www.sbml.org/sbml/level3/version1/comp/version1" level="3" version="1" comp:required="true">
<model id="m"><comp:listOfSubmodels><comp:submodel comp:id="a" comp:modelRef="b"/></comp:listOfSubmodels></model></sbml>"#;
    assert!(KineticModel::from_sbml_str(text).unwrap_err().to_string().contains("comp"));
    // Declaring it without using it changes nothing.
    let unused = text
        .replace(r#"<comp:listOfSubmodels><comp:submodel comp:id="a" comp:modelRef="b"/></comp:listOfSubmodels>"#, "");
    assert!(KineticModel::from_sbml_str(&unused).is_ok());
}

fn try_ensemble(name: &str, end: f64, runs: usize, algorithm: Algorithm) -> Result<f64, String> {
    let m = model(name);
    let e = ensemble(&m, &[0.0, end], runs, algorithm, 3).map_err(|e| e.to_string())?;
    let symbol = if m.index.contains_key("X") { "X" } else { "S" };
    Ok(e.mean[1][m.slot(symbol).unwrap()])
}

#[test]
fn leaps_are_bounded_by_autocatalysts_and_catalysts() {
    // X -> 2X with k = 1 from 1000: mean 1000 e^5. Explicit leaping is low by about 7 % at epsilon 0.03.
    let exact = 1000.0 * 5f64.exp();
    let m = model("yule");
    for steps in [1, 10, 100] {
        let (mean, _) = mean_at_end(&m, "X", 5.0, steps, 100, TAU);
        assert!((mean - exact).abs() < 0.12 * exact, "{steps} stops: {mean} vs {exact}");
    }
    let (mean, _) = mean_at_end(&m, "X", 5.0, 1, 40, Algorithm::TauLeap { epsilon: 0.001 });
    assert!((mean - exact).abs() < 0.03 * exact, "{mean}");
    // A + B -> 2B is logistic: B(5) is 12993.
    let (mean, _) = mean_at_end(&model("si"), "B", 5.0, 1, 100, TAU);
    assert!((mean - 12993.0).abs() < 0.1 * 12993.0, "{mean}");
    // A + E -> B + E, with E decaying: A(1) of about 27434 by the ODE.
    let (mean, _) = mean_at_end(&model("cat"), "A", 1.0, 1, 40, TAU);
    assert!((mean - 27434.0).abs() < 0.05 * 27434.0, "{mean}");
}

#[test]
fn a_rate_that_turns_on_just_after_a_stop_is_read_there() {
    for name in ["win_gtleq", "win_geqlt"] {
        for algorithm in [Algorithm::Direct, TAU] {
            let mean = try_ensemble(name, 1000.0, 200, algorithm).unwrap();
            assert!((mean - 100.0).abs() < 5.0, "{name} {algorithm:?}: {mean}");
        }
    }
}

#[test]
fn stochastic_pulses_from_rounding_are_counted_and_smooth_forcing_is_refused() {
    // 100 per unit time for one unit in every ten.
    let mean = try_ensemble("c_pul", 1000.0, 20, Algorithm::Direct).unwrap();
    assert!((mean - 10_000.0).abs() < 300.0, "{mean}");
    let error = try_ensemble("c_sin", 100.0, 2, Algorithm::Direct).unwrap_err();
    assert!(error.contains("continuously"), "{error}");
}

#[test]
fn stochastic_simulation_counts_whole_molecules() {
    assert!(try_ensemble("frac", 1.0, 2, Algorithm::Direct).unwrap_err().contains("whole number"));
    assert!(try_ensemble("frac2", 1.0, 2, Algorithm::Direct).unwrap_err().contains("whole number"));
}

#[test]
fn an_explosive_model_stops_with_a_reason() {
    let error = try_ensemble("ssa_expl", 1.0, 1, TAU).unwrap_err();
    assert!(error.contains("without limit"), "{error}");
}

#[test]
fn newton_does_not_accept_a_root_at_infinity() {
    for name in ["ss_inf", "ss_inf2"] {
        assert!(steady_state(&model(name), &SteadyOptions::default()).is_err(), "{name}");
    }
}

#[test]
fn steady_states_look_past_events_and_hidden_jumps() {
    let value = |name: &str| {
        let m = model(name);
        let s = steady_state(&m, &SteadyOptions::default()).unwrap();
        s.values[m.slot("S").unwrap()]
    };
    assert!((value("ss_ev") - 5.0).abs() < 1e-9, "event after the horizon");
    assert!((value("o_hid") - 5.0).abs() < 1e-9, "time - 1000 > 0");
    let relaxation = steady_state(&model("ss_relax"), &SteadyOptions::default()).unwrap_err().to_string();
    assert!(relaxation.contains("fires"), "{relaxation}");
    let periodic = steady_state(&model("o_per"), &SteadyOptions::default()).unwrap_err().to_string();
    assert!(periodic.contains("depend on the time"), "{periodic}");
}

#[test]
fn a_slow_instability_is_reported_and_an_empty_model_is_stable() {
    assert!(!steady_state(&model("ss_slow_unst"), &SteadyOptions::default()).unwrap().stable);
    for name in ["ss_empty", "ss_empty2"] {
        let m = model(name);
        let s = steady_state(&m, &SteadyOptions::default()).unwrap();
        assert!(s.stable, "{name}");
        control_analysis(&m, &s).unwrap();
    }
}

#[test]
fn a_delay_decides_its_own_stability() {
    // y' = -y(t - tau) is stable up to tau = pi / 2.
    let stable = steady_state(&model("f_dde_1.0"), &SteadyOptions::default()).unwrap();
    assert!(stable.stable);
    let unstable = steady_state(&model("f_dde_2.0"), &SteadyOptions::default()).unwrap();
    assert!(!unstable.stable);
    assert!(unstable.largest_eigenvalue.is_nan());
}

#[test]
fn local_sensitivities_see_parameters_that_set_the_start() {
    // C = Etot kon S / (kon S + koff + kcat), so dC/dEtot = 5/6.
    let m = model("sens_moiety");
    let state = steady_state(&m, &SteadyOptions::default()).unwrap();
    let control = control_analysis(&m, &state).unwrap();
    let local = local_steady(&m, &state, &control, &["Etot".to_string()]).unwrap();
    let c = control.states.iter().position(|s| s == "C#a").unwrap();
    assert!((local.state[(c, 0)] - 5.0 / 6.0).abs() < 1e-7, "{}", local.state[(c, 0)]);
    let local = sensitivities("sens_bnd", &["p"]);
    assert!((local.flux[(0, 0)] - 2.0).abs() < 1e-7, "{}", local.flux[(0, 0)]);
    let local = sensitivities("sens_comp", &["c"]);
    assert!((local.flux[(0, 0)] - 1.0).abs() < 1e-7, "{}", local.flux[(0, 0)]);
}

fn profile_of(
    name: &str,
    sigma: f64,
    start: f64,
    values: &[(f64, f64)],
    symbol: &str,
    factor: f64,
    points: usize,
) -> (rashevsky::kinetic::fit::Profile, f64) {
    let m = model(name);
    let data: Vec<Observation> =
        values.iter().map(|&(time, value)| Observation { time, symbol: symbol.into(), value, sigma }).collect();
    let parameters = [FitParameter { name: "k".into(), start }];
    let options = SimOptions::default();
    let fitted = fit(&m, &parameters, &data, &options).unwrap();
    (profile(&m, &parameters, &data, &fitted, "k", factor, points, &options).unwrap(), fitted.estimates[0])
}

#[test]
fn profile_limits_stay_right_when_the_likelihood_is_sharp() {
    // A' = -k A from 2. chi2(k) has a closed form here, so the limits are checked against its roots.
    let times = [0.5, 1.0, 2.0, 3.0];
    let curve = |k: f64, t: f64| 2.0 * (-k * t).exp();
    let observed: Vec<(f64, f64)> = times.iter().map(|&t| (t, curve(0.7, t))).collect();
    for sigma in [0.02, 0.002, 0.0002] {
        let (p, _) = profile_of("p_decay", sigma, 1.0, &observed, "A", 4.0, 15);
        let chi2 = |k: f64| -> f64 { times.iter().map(|&t| ((curve(k, t) - curve(0.7, t)) / sigma).powi(2)).sum() };
        for limit in [p.lower.unwrap(), p.upper.unwrap()] {
            assert!((chi2(limit) - 3.841_458_820_694_124).abs() < 0.02, "sigma {sigma}: {limit} gives {}", chi2(limit));
        }
    }
}

#[test]
fn a_profile_reports_a_fit_stuck_in_a_local_minimum() {
    // S' = cos(k t) with data from k = 1, and a fit started beyond the first minimum.
    let times: Vec<f64> = (1..=12).map(|i| f64::from(i) * 0.5).collect();
    let observed: Vec<(f64, f64)> = times.iter().map(|&t| (t, t.sin())).collect();
    let (p, estimate) = profile_of("p_osc", 0.05, 3.0, &observed, "S", 8.0, 31);
    assert!((estimate - 1.0).abs() > 0.5, "this start should stick: {estimate}");
    assert!((p.lowest_at - 1.0).abs() < 0.2, "{}", p.lowest_at);
    assert!(p.lowest_chi2 < 1000.0, "{}", p.lowest_chi2);
}

#[test]
fn every_way_of_writing_a_time_condition_is_found() {
    for method in [Method::Rodas4, Method::Bdf] {
        for name in ["th_narylt50time5001", "th_2time100", "th_abstimecw", "th_rule", "th_fn", "th_par"] {
            let s = end_value(&model(name), "S", 100.0, 1, method);
            assert!((s - 1.0).abs() < 1e-5, "{name} {method:?}: {s}");
        }
        let s = end_value(&model("th_per"), "S", 100.0, 1, method);
        assert!((s - 10.0).abs() < 1e-5, "periodic {method:?}: {s}");
        for name in ["tev_nary", "tev_geq"] {
            let f = end_value(&model(name), "f", 100.0, 1, method);
            assert_eq!(f, 1.0, "{name} {method:?}");
        }
    }
}

#[test]
fn a_late_switch_does_not_stop_either_integrator() {
    for method in [Method::Rodas4, Method::Bdf] {
        let s = end_value(&model("o_th"), "S", 2e9, 1, method);
        assert!((s / 1e9 - 1.0).abs() < 1e-6, "{method:?}: {s}");
    }
}

#[test]
fn overrides_refuse_what_rules_set_and_replace_their_alias() {
    let mut m = model("set_assigned");
    assert!(m.set_value("S#a", 100.0).is_err());
    assert!(m.set_value("S", 100.0).is_err());
    let mut m = model("set_alg");
    assert!(m.set_value("S#a", 100.0).is_err());
    let mut a = model("set_two");
    a.set_value("A", 5.0).unwrap();
    a.set_value("A#a", 7.0).unwrap();
    assert_eq!(a.initial[a.slot("A#a").unwrap()], 7.0);
    let mut b = model("set_two");
    b.set_value("A#a", 7.0).unwrap();
    b.set_value("A", 5.0).unwrap();
    assert_eq!(b.initial[b.slot("A").unwrap()], 5.0);
}

#[test]
fn a_fast_scale_follows_overrides() {
    let mut m = model("fast_1");
    m.set_value("kf", 1e-9).unwrap();
    m.set_value("kr", 1e-9).unwrap();
    let run = simulate(&m, &[0.0], &SimOptions::default()).unwrap();
    assert!((run.series(&m, "A").unwrap()[0] - 5.0).abs() < 1e-3, "{:?}", run.series(&m, "A").unwrap());
}

#[test]
fn sobol_of_an_unresponsive_output_is_zero_and_unstable_samples_do_not_count() {
    use rashevsky::kinetic::sensitivity::{Output, Range, sobol_model};
    let m = model("o_bruss");
    let output = Output::Steady { symbol: "X".into(), integrate: 100.0 };
    // Below b = 2 the fixed point attracts and X = a whatever b is.
    let quiet = [Range { name: "b".into(), low: 1.0, high: 1.5, log: false }];
    let r = sobol_model(&m, &quiet, 16, 1, &output).unwrap();
    assert_eq!(r.first, vec![0.0]);
    assert_eq!(r.total, vec![0.0]);
    // Above it there is no attractor to read.
    let wide = [Range { name: "b".into(), low: 1.0, high: 3.0, log: false }];
    assert!(sobol_model(&m, &wide, 16, 1, &output).is_err());
}
