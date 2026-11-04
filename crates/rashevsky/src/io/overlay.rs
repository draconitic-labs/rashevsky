//! Flux overlays on the metabolism map of the Biochemistry site (`biochemistry.draconitic.com`).
//!
//! The map draws Reactome events. A model reaction reaches an event through a shared Reactome
//! identifier, through a Rhea identifier the event carries, or through a shared EC number. Each overlay
//! entry says which: `reactome` and `rhea` are identifier matches, while `ec` is an approximation,
//! since one EC number can cover several events and compartments. A model and the map are separate
//! reconstructions, so the overlay is a projection for reading fluxes in context; it does not assert
//! that an event carries the given flux in human cells.

use crate::error::Result;
use crate::model::Model;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;

/// Bump when the overlay's fields change meaning.
pub const OVERLAY_FORMAT: &str = "rashevsky-flux-overlay/1";

/// A node of the scene.
#[derive(Clone, Debug, Deserialize)]
pub struct SceneNode {
    /// Horizontal position.
    pub x: f64,
    /// Vertical position.
    pub y: f64,
    /// Short label.
    #[serde(default)]
    pub l: Option<String>,
    /// Hub flag.
    #[serde(default)]
    pub h: f64,
}

/// A drawn reaction of the scene.
#[derive(Clone, Debug, Deserialize)]
pub struct SceneReaction {
    /// Horizontal position.
    pub x: f64,
    /// Vertical position.
    pub y: f64,
    /// Input nodes.
    pub i: Vec<usize>,
    /// Output nodes.
    pub o: Vec<usize>,
    /// Reactome events drawn here.
    pub e: Vec<String>,
    /// EC numbers.
    #[serde(default)]
    pub ec: Vec<String>,
}

/// A process of the scene, drawn as a point.
#[derive(Clone, Debug, Deserialize)]
pub struct SceneProcess {
    /// Horizontal position.
    pub x: f64,
    /// Vertical position.
    pub y: f64,
    /// Reactome events.
    pub e: Vec<String>,
    /// EC numbers.
    #[serde(default)]
    pub ec: Vec<String>,
}

/// World size of the scene.
#[derive(Clone, Copy, Debug, Deserialize)]
pub struct World {
    /// Width.
    pub w: f64,
    /// Height.
    pub h: f64,
}

/// A domain label of the scene.
#[derive(Clone, Debug, Deserialize)]
pub struct SceneDomain {
    /// Label.
    pub label: String,
    /// Colour.
    pub color: String,
    /// Horizontal position.
    pub x: f64,
    /// Vertical position.
    pub y: f64,
}

/// The fields of the site's scene file that the overlay needs.
#[derive(Clone, Debug, Deserialize)]
pub struct Scene {
    /// Release identifier.
    pub release: String,
    /// World size.
    pub world: World,
    /// Domains.
    pub domains: Vec<SceneDomain>,
    /// Metabolite nodes.
    pub nodes: Vec<SceneNode>,
    /// Drawn reactions.
    pub reactions: Vec<SceneReaction>,
    /// Processes.
    pub processes: Vec<SceneProcess>,
}

/// Flux assigned to one event.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct OverlayEvent {
    /// Sum of absolute fluxes of the matched reactions, mmol/gDW/h.
    pub flux: f64,
    /// Strongest match type: `reactome`, `rhea` or `ec`.
    #[serde(rename = "match")]
    pub match_type: String,
    /// Model reactions that contribute.
    pub reactions: Vec<String>,
}

/// Coverage of an overlay.
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize)]
pub struct Coverage {
    /// Events in the scene.
    pub events_in_scene: usize,
    /// Events reached by a Reactome or Rhea identifier.
    pub events_exact: usize,
    /// Events reached only by an EC number.
    pub events_ec_only: usize,
    /// Events that carry non-zero flux.
    pub events_with_flux: usize,
    /// Model reactions matched to at least one event.
    pub reactions_matched: usize,
    /// Model reactions that carry flux.
    pub reactions_with_flux: usize,
}

/// A flux overlay.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Overlay {
    /// Format identifier.
    pub format: String,
    /// Scene release the event identifiers belong to.
    pub scene_release: String,
    /// Model identifier.
    pub model: String,
    /// Objective value of the flux solution, when given.
    pub objective: Option<f64>,
    /// Flux per Reactome event.
    pub events: BTreeMap<String, OverlayEvent>,
    /// Coverage counts.
    pub coverage: Coverage,
}

fn annotation_value<'a>(resource: &'a str, database: &str) -> Option<&'a str> {
    let marker = format!("identifiers.org/{database}/");
    if let Some(at) = resource.find(&marker) {
        return Some(&resource[at + marker.len()..]);
    }
    let colon = format!("identifiers.org/{database}:");
    resource.find(&colon).map(|at| &resource[at + colon.len()..])
}

/// Projects a flux vector onto the scene's events.
pub fn build_overlay(
    model: &Model,
    fluxes: &[f64],
    objective: Option<f64>,
    scene: &Scene,
    rhea_events: &BTreeMap<String, Vec<String>>,
) -> Result<Overlay> {
    let mut known: HashMap<&str, ()> = HashMap::new();
    let mut by_ec: HashMap<&str, Vec<&str>> = HashMap::new();
    for (events, ecs) in
        scene.reactions.iter().map(|r| (&r.e, &r.ec)).chain(scene.processes.iter().map(|p| (&p.e, &p.ec)))
    {
        for e in events {
            known.insert(e.as_str(), ());
        }
        for ec in ecs {
            by_ec.entry(ec.as_str()).or_default().extend(events.iter().map(String::as_str));
        }
    }
    // Per event, reactions reached by an identifier and reactions reached only through an EC number.
    // An event with any identifier match reports those reactions alone, so that an approximate
    // association never adds to, or hides inside, an exact one.
    #[derive(Default)]
    struct Tally {
        reactome: Vec<usize>,
        rhea: Vec<usize>,
        ec: Vec<usize>,
    }
    let mut tally: BTreeMap<String, Tally> = BTreeMap::new();
    let mut reactions_matched = 0;
    let mut reactions_with_flux = 0;
    for (j, r) in model.reactions.iter().enumerate() {
        if fluxes[j].abs() > 1e-9 {
            reactions_with_flux += 1;
        }
        let mut exact = false;
        let mut any = false;
        for resource in &r.annotations {
            if let Some(id) = annotation_value(resource, "reactome")
                && known.contains_key(id)
            {
                tally.entry(id.to_string()).or_default().reactome.push(j);
                exact = true;
            }
            if let Some(id) = annotation_value(resource, "rhea")
                && let Some(list) = rhea_events.get(id)
            {
                for e in list.iter().filter(|e| known.contains_key(e.as_str())) {
                    tally.entry(e.clone()).or_default().rhea.push(j);
                    exact = true;
                }
            }
        }
        any |= exact;
        if !exact {
            for resource in &r.annotations {
                if let Some(ec) = annotation_value(resource, "ec-code")
                    && let Some(list) = by_ec.get(ec)
                {
                    for e in list {
                        tally.entry((*e).to_string()).or_default().ec.push(j);
                        any = true;
                    }
                }
            }
        }
        if any {
            reactions_matched += 1;
        }
    }
    let mut events: BTreeMap<String, OverlayEvent> = BTreeMap::new();
    for (event, t) in tally {
        let (kind, mut members) = if !t.reactome.is_empty() {
            let mut m = t.reactome;
            m.extend(t.rhea);
            ("reactome", m)
        } else if !t.rhea.is_empty() {
            ("rhea", t.rhea)
        } else {
            ("ec", t.ec)
        };
        members.sort_unstable();
        members.dedup();
        events.insert(
            event,
            OverlayEvent {
                flux: members.iter().map(|&j| fluxes[j].abs()).sum(),
                match_type: kind.to_string(),
                reactions: members.iter().map(|&j| model.reactions[j].id.clone()).collect(),
            },
        );
    }
    let coverage = Coverage {
        events_in_scene: known.len(),
        events_exact: events.values().filter(|e| e.match_type != "ec").count(),
        events_ec_only: events.values().filter(|e| e.match_type == "ec").count(),
        events_with_flux: events.values().filter(|e| e.flux > 1e-9).count(),
        reactions_matched,
        reactions_with_flux,
    };
    Ok(Overlay {
        format: OVERLAY_FORMAT.to_string(),
        scene_release: scene.release.clone(),
        model: model.id.clone(),
        objective,
        events,
        coverage,
    })
}

/// Colour on a logarithmic flux scale, from teal through amber to rose.
fn flux_color(t: f64) -> String {
    let stops: [(f64, [f64; 3]); 3] =
        [(0.0, [46.0, 196.0, 182.0]), (0.5, [244.0, 185.0, 66.0]), (1.0, [240.0, 84.0, 110.0])];
    let t = t.clamp(0.0, 1.0);
    let (a, b) = if t < 0.5 { (stops[0], stops[1]) } else { (stops[1], stops[2]) };
    let w = (t - a.0) / (b.0 - a.0);
    let c = |k: usize| (a.1[k] + w * (b.1[k] - a.1[k])).round() as u8;
    format!("#{:02x}{:02x}{:02x}", c(0), c(1), c(2))
}

/// Draws the overlay on the scene as a standalone SVG.
pub fn render_svg(scene: &Scene, overlay: &Overlay) -> String {
    let max_flux = overlay.events.values().map(|e| e.flux).fold(0.0f64, f64::max).max(1e-9);
    let scale = |flux: f64| -> f64 {
        ((flux.max(1e-3) / 1e-3).ln() / (max_flux.max(1e-3) / 1e-3).ln().max(1e-9)).clamp(0.0, 1.0)
    };
    let (w, h) = (scene.world.w, scene.world.h);
    let mut s = String::with_capacity(1 << 20);
    let _ = write!(
        s,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {w} {h}\" font-family=\"Inter,Helvetica,Arial,sans-serif\">\
         <rect width=\"{w}\" height=\"{h}\" fill=\"#0e1013\"/>"
    );
    for d in &scene.domains {
        let label = crate::sbml::write::escape(&d.label);
        let _ = write!(
            s,
            "<text x=\"{}\" y=\"{}\" fill=\"{}\" font-size=\"64\" opacity=\"0.5\" text-anchor=\"middle\">{label}</text>",
            d.x, d.y, d.color
        );
    }
    // Reactions without flux first, so that flux draws on top.
    let event_flux = |events: &[String]| -> Option<(f64, &str)> {
        let mut total = 0.0;
        let mut best: Option<&str> = None;
        for e in events {
            if let Some(o) = overlay.events.get(e) {
                total += o.flux;
                if best.is_none() || o.match_type != "ec" {
                    best = Some(o.match_type.as_str());
                }
            }
        }
        best.filter(|_| total > 1e-9).map(|m| (total, m))
    };
    for pass in [false, true] {
        for r in &scene.reactions {
            let flux = event_flux(&r.e);
            if flux.is_some() != pass {
                continue;
            }
            let (stroke, width, opacity, dash) = match flux {
                Some((f, m)) => {
                    let t = scale(f);
                    (flux_color(t), 4.0 + 14.0 * t, 0.95, if m == "ec" { " stroke-dasharray=\"12 8\"" } else { "" })
                }
                None => ("#3a404a".to_string(), 2.0, 0.7, ""),
            };
            for &k in r.i.iter().chain(&r.o) {
                if let Some(n) = scene.nodes.get(k) {
                    let _ = write!(
                        s,
                        "<line x1=\"{:.0}\" y1=\"{:.0}\" x2=\"{:.0}\" y2=\"{:.0}\" stroke=\"{stroke}\" stroke-width=\"{width:.1}\" opacity=\"{opacity}\"{dash}/>",
                        r.x, r.y, n.x, n.y
                    );
                }
            }
        }
    }
    for n in &scene.nodes {
        let _ = write!(
            s,
            "<circle cx=\"{:.0}\" cy=\"{:.0}\" r=\"{}\" fill=\"#8a93a2\"/>",
            n.x,
            n.y,
            if n.h > 0.0 { 9 } else { 5 }
        );
    }
    for p in &scene.processes {
        if let Some((f, _)) = event_flux(&p.e) {
            let _ =
                write!(s, "<circle cx=\"{:.0}\" cy=\"{:.0}\" r=\"9\" fill=\"{}\"/>", p.x, p.y, flux_color(scale(f)));
        }
    }
    s.push_str("</svg>\n");
    s
}
