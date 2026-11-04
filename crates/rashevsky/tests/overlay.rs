//! How model reactions reach the events of the metabolism map, and what the drawing says about them.

use rashevsky::Model;
use rashevsky::io::overlay::{Scene, build_overlay, render_svg};
use rashevsky::model::Reaction;
use std::collections::BTreeMap;

/// One drawn reaction over the given events, which all carry the EC number 1.1.1.1.
fn scene(events: &[&str]) -> Scene {
    serde_json::from_value(serde_json::json!({
        "release": "test",
        "world": {"w": 100.0, "h": 100.0},
        "domains": [],
        "nodes": [{"x": 0.0, "y": 0.0}, {"x": 10.0, "y": 0.0}],
        "reactions": [{"x": 5.0, "y": 0.0, "i": [0], "o": [1], "e": events, "ec": ["1.1.1.1"]}],
        "processes": [],
    }))
    .unwrap()
}

fn model(annotations: &[&[&str]]) -> Model {
    let mut m = Model::new("overlay");
    for (k, list) in annotations.iter().enumerate() {
        m.reactions.push(Reaction {
            id: format!("r{k}"),
            annotations: list.iter().map(|s| (*s).to_string()).collect(),
            ..Reaction::default()
        });
    }
    m.reindex();
    m
}

const NO_RHEA: BTreeMap<String, Vec<String>> = BTreeMap::new();

#[test]
fn both_reactome_namespaces_match() {
    let m = model(&[
        &["http://identifiers.org/reactome/R-HSA-1"],
        &["http://identifiers.org/reactome.reaction/R-HSA-2"],
        &["https://identifiers.org/reactome.reaction:R-HSA-3"],
    ]);
    let s = scene(&["R-HSA-1", "R-HSA-2", "R-HSA-3"]);
    let overlay = build_overlay(&m, &[1.0, 2.0, 3.0], None, &s, &NO_RHEA).unwrap();
    assert_eq!(overlay.events["R-HSA-1"].flux, 1.0);
    assert_eq!(overlay.events["R-HSA-2"].flux, 2.0);
    assert_eq!(overlay.events["R-HSA-3"].flux, 3.0);
    assert_eq!(overlay.events["R-HSA-2"].match_type, "reactome");
}

#[test]
fn an_exact_match_without_flux_does_not_hide_an_ec_match_with_flux() {
    let m = model(&[&["http://identifiers.org/reactome/R-HSA-1"], &["http://identifiers.org/ec-code/1.1.1.1"]]);
    let s = scene(&["R-HSA-1"]);
    let overlay = build_overlay(&m, &[0.0, 5.0], None, &s, &NO_RHEA).unwrap();
    assert_eq!(overlay.events["R-HSA-1"].match_type, "ec");
    assert_eq!(overlay.events["R-HSA-1"].flux, 5.0);
    assert_eq!(overlay.coverage.events_with_flux, 1);
}

#[test]
fn an_exact_match_with_flux_reports_alone() {
    let m = model(&[&["http://identifiers.org/reactome/R-HSA-1"], &["http://identifiers.org/ec-code/1.1.1.1"]]);
    let s = scene(&["R-HSA-1"]);
    let overlay = build_overlay(&m, &[2.0, 5.0], None, &s, &NO_RHEA).unwrap();
    assert_eq!(overlay.events["R-HSA-1"].match_type, "reactome");
    assert_eq!(overlay.events["R-HSA-1"].flux, 2.0);
}

#[test]
fn one_reaction_on_several_events_is_drawn_at_its_own_value() {
    let ec = model(&[&["http://identifiers.org/ec-code/1.1.1.1"]]);
    let four = scene(&["R-HSA-1", "R-HSA-2", "R-HSA-3", "R-HSA-4"]);
    let one = scene(&["R-HSA-1"]);
    let drawn_four = render_svg(&four, &build_overlay(&ec, &[5.0], None, &four, &NO_RHEA).unwrap());
    let drawn_one = render_svg(&one, &build_overlay(&ec, &[5.0], None, &one, &NO_RHEA).unwrap());
    assert_eq!(drawn_four, drawn_one);
}

#[test]
fn approximate_flux_is_drawn_as_approximate() {
    // The exact event carries nothing, so the flux shown comes from the EC match and must look it.
    let m = model(&[&["http://identifiers.org/reactome/R-HSA-1"], &["http://identifiers.org/ec-code/1.1.1.1"]]);
    let s = scene(&["R-HSA-1", "R-HSA-2"]);
    let approximate = render_svg(&s, &build_overlay(&m, &[0.0, 5.0], None, &s, &NO_RHEA).unwrap());
    assert!(approximate.contains("stroke-dasharray"), "{approximate}");
    let exact = render_svg(&s, &build_overlay(&m, &[5.0, 0.0], None, &s, &NO_RHEA).unwrap());
    assert!(!exact.contains("stroke-dasharray"), "{exact}");
}
