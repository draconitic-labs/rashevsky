//! Checks of the kinetic engine on analytic cases and on the bundled glycolysis model.

use rashevsky::kinetic::ssa::{Algorithm, ensemble};
use rashevsky::kinetic::{KineticModel, Method, SimOptions, SteadyOptions, control_analysis, simulate, steady_state};

const TEUSINK: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../models/teusink/BIOMD0000000064.xml");

fn teusink() -> KineticModel {
    KineticModel::read_sbml(TEUSINK).unwrap()
}

/// `A -> B` with mass action, and `B <-> C` near equilibrium.
const CHAIN: &str = r#"<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
<model id="chain"><listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
<listOfSpecies>
<species id="A" compartment="c" initialConcentration="2" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
<species id="B" compartment="c" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
</listOfSpecies>
<listOfParameters><parameter id="k" value="0.7" constant="true"/></listOfParameters>
<listOfReactions><reaction id="r" reversible="false">
<listOfReactants><speciesReference species="A" stoichiometry="1" constant="true"/></listOfReactants>
<listOfProducts><speciesReference species="B" stoichiometry="1" constant="true"/></listOfProducts>
<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci>k</ci><ci>A</ci><ci>c</ci></apply></math></kineticLaw>
</reaction></listOfReactions></model></sbml>"#;

#[test]
fn analytic_decay_at_both_integrators() {
    let m = KineticModel::from_sbml_str(CHAIN).unwrap();
    let times = [0.0, 0.5, 1.0, 2.0, 4.0];
    for method in [Method::Rodas4, Method::Bdf] {
        let options = SimOptions { method, rtol: 1e-10, atol: 1e-13, ..SimOptions::default() };
        let run = simulate(&m, &times, &options).unwrap();
        let a = run.series(&m, "A").unwrap();
        let b = run.series(&m, "B").unwrap();
        for (k, &t) in times.iter().enumerate() {
            let exact = 2.0 * (-0.7 * t).exp();
            assert!((a[k] - exact).abs() < 2e-8, "{method:?} t={t}: {} vs {exact}", a[k]);
            // Conservation of A + B.
            assert!((a[k] + b[k] - 2.0).abs() < 1e-10);
        }
    }
}

#[test]
fn analytic_jacobian_matches_finite_differences() {
    let m = teusink();
    let run = simulate(&m, &[0.0, 5.0], &SimOptions::default()).unwrap();
    let values = &run.values[1];
    let y: Vec<f64> = m.state_slots.iter().map(|&s| values[s]).collect();
    let n = y.len();
    use rashevsky::kinetic::ode::OdeSystem;
    use rashevsky::kinetic::sim::KineticSystem;
    let analytic = KineticSystem::new(&m, false);
    let numeric = KineticSystem::new(&m, true);
    let mut ja = rashevsky::nalgebra::DMatrix::zeros(n, n);
    let mut jn = rashevsky::nalgebra::DMatrix::zeros(n, n);
    analytic.jacobian(5.0, &y, &mut ja);
    numeric.jacobian(5.0, &y, &mut jn);
    let scale = ja.iter().fold(0.0f64, |a, v| a.max(v.abs()));
    for i in 0..n {
        for j in 0..n {
            assert!(
                (ja[(i, j)] - jn[(i, j)]).abs() < 2e-4 * scale.max(1.0),
                "J[{i},{j}]: {} vs {}",
                ja[(i, j)],
                jn[(i, j)]
            );
        }
    }
}

#[test]
fn steady_state_has_zero_rates_and_keeps_moieties() {
    let m = teusink();
    let s = steady_state(&m, &SteadyOptions::default()).unwrap();
    assert!(s.residual < 1e-9, "{}", s.residual);
    let value = |name: &str| s.values[m.slot(name).unwrap()];
    // NAD + NADH is conserved at its initial total, 1.59.
    assert!((value("NAD") + value("NADH") - 1.59).abs() < 1e-8);
}

#[test]
fn control_theorems_hold_on_glycolysis() {
    let m = teusink();
    let s = steady_state(&m, &SteadyOptions::default()).unwrap();
    let c = control_analysis(&m, &s).unwrap();
    let t = c.theorems;
    for (name, value) in [
        ("flux summation", t.flux_summation),
        ("concentration summation", t.concentration_summation),
        ("flux connectivity", t.flux_connectivity),
        ("concentration connectivity", t.concentration_connectivity),
        ("flux kernel", t.flux_kernel),
        ("concentration kernel", t.concentration_kernel),
    ] {
        assert!(value < 1e-8, "{name}: {value}");
    }
}

#[test]
fn control_coefficient_matches_a_perturbation() {
    // C^J_i = d ln J / d ln e_i, tested by scaling the transporter's rate through its Vmax.
    let m = teusink();
    let s = steady_state(&m, &SteadyOptions::default()).unwrap();
    let c = control_analysis(&m, &s).unwrap();
    let glt = c.reactions.iter().position(|r| r == "vGLT").unwrap();
    let adh = c.reactions.iter().position(|r| r == "vADH").unwrap();
    let analytic = c.scaled_flux_control[(adh, glt)];
    // Scale the transport law directly by compiling a copy with a modified Vmax.
    let text = std::fs::read_to_string(TEUSINK).unwrap();
    let bump = 1.001f64;
    let needle = "id=\"VmGLT\" value=\"97.264\"";
    assert!(text.contains(needle), "the local Vmax of the transporter is expected in this form");
    let high = KineticModel::from_sbml_str(&text.replace(needle, &format!("id=\"VmGLT\" value=\"{}\"", 97.264 * bump)))
        .unwrap();
    let low = KineticModel::from_sbml_str(&text.replace(needle, &format!("id=\"VmGLT\" value=\"{}\"", 97.264 / bump)))
        .unwrap();
    let flux = |model: &KineticModel| steady_state(model, &SteadyOptions::default()).unwrap().fluxes[adh];
    let numeric = (flux(&high).ln() - flux(&low).ln()) / (2.0 * bump.ln());
    assert!((analytic - numeric).abs() < 1e-3, "{analytic} vs {numeric}");
}

#[test]
fn stochastic_means_match_the_deterministic_decay() {
    let text = CHAIN
        .replace("initialConcentration=\"2\"", "initialAmount=\"200\"")
        .replace("initialConcentration=\"0\"", "initialAmount=\"0\"");
    let m = KineticModel::from_sbml_str(&text).unwrap();
    let times = [0.0, 1.0, 2.0];
    let exact = |t: f64| 200.0 * (-0.7 * t).exp();
    let slot = m.slot("A#a").unwrap();
    for algorithm in [Algorithm::Direct, Algorithm::TauLeap { epsilon: 0.03 }] {
        let e = ensemble(&m, &times, 4000, algorithm, 3).unwrap();
        for (k, &t) in times.iter().enumerate() {
            // Binomial variance of pure death.
            let p = (-0.7 * t).exp();
            let se = (200.0 * p * (1.0 - p) / 4000.0).sqrt().max(1e-6);
            assert!((e.mean[k][slot] - exact(t)).abs() < 5.0 * se + 0.05, "{algorithm:?} t={t}: {}", e.mean[k][slot]);
        }
    }
}

#[test]
fn fit_recovers_a_rate_and_profile_bounds_it() {
    use rashevsky::kinetic::fit::{FitParameter, Observation, fit, profile};
    let m = KineticModel::from_sbml_str(CHAIN).unwrap();
    // Noise-free data with a stated measurement error, to test the machinery.
    let data: Vec<Observation> = [0.5, 1.0, 2.0, 3.0]
        .iter()
        .map(|&t| Observation { time: t, symbol: "A".into(), value: 2.0 * (-0.7f64 * t).exp(), sigma: 0.02 })
        .collect();
    let parameters = [FitParameter { name: "k".into(), start: 3.0 }];
    let options = SimOptions::default();
    let r = fit(&m, &parameters, &data, &options).unwrap();
    assert!((r.estimates[0] - 0.7).abs() < 1e-4, "{:?}", r.estimates);
    assert!(r.chi2 < 1e-6);
    let p = profile(&m, &parameters, &data, &r, "k", 1.5, 21, &options).unwrap();
    let (low, high) = (p.lower.unwrap(), p.upper.unwrap());
    assert!(low < 0.7 && high > 0.7, "{low} {high}");
    // With sigma 0.02 the interval is narrow.
    assert!(high / low < 1.2, "{low} {high}");
}

#[test]
fn a_non_identifiable_pair_has_an_open_profile() {
    use rashevsky::kinetic::fit::{FitParameter, Observation, fit, profile};
    // Only the product of a and b enters the rate, so neither is identifiable alone.
    let text = CHAIN
        .replace(
            "<parameter id=\"k\" value=\"0.7\" constant=\"true\"/>",
            "<parameter id=\"a\" value=\"1\" constant=\"true\"/><parameter id=\"b\" value=\"0.7\" constant=\"true\"/>",
        )
        .replace("<ci>k</ci>", "<ci>a</ci><ci>b</ci>");
    let m = KineticModel::from_sbml_str(&text).unwrap();
    let data: Vec<Observation> = [0.5, 1.0, 2.0, 3.0]
        .iter()
        .map(|&t| Observation { time: t, symbol: "A".into(), value: 2.0 * (-0.7f64 * t).exp(), sigma: 0.02 })
        .collect();
    let parameters = [FitParameter { name: "a".into(), start: 1.5 }, FitParameter { name: "b".into(), start: 1.0 }];
    let options = SimOptions::default();
    let r = fit(&m, &parameters, &data, &options).unwrap();
    assert!((r.estimates[0] * r.estimates[1] - 0.7).abs() < 1e-3, "{:?}", r.estimates);
    let p = profile(&m, &parameters, &data, &r, "a", 3.0, 9, &options).unwrap();
    assert!(p.lower.is_none() && p.upper.is_none(), "{p:?}");
}

#[test]
fn a_value_override_changes_the_simulation() {
    let mut m = KineticModel::from_sbml_str(CHAIN).unwrap();
    m.set_value("k", 1.4).unwrap();
    let run = simulate(&m, &[0.0, 1.0], &SimOptions::default()).unwrap();
    let a = run.series(&m, "A").unwrap();
    assert!((a[1] - 2.0 * (-1.4f64).exp()).abs() < 1e-7);
}
