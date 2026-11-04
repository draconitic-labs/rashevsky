//! SBML document reader.

use super::dom::{self, Node};
use super::mathml::{read_lambda, read_math};
use super::*;
use crate::error::{Error, Result};
use std::io::Read as _;
use std::path::Path;

/// Reads an SBML file; `.gz` files are decompressed.
pub fn read_file(path: impl AsRef<Path>) -> Result<Document> {
    let path = path.as_ref();
    let bytes = std::fs::read(path)?;
    let text = if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut out = String::new();
        flate2::read::GzDecoder::new(&bytes[..]).read_to_string(&mut out)?;
        out
    } else {
        String::from_utf8(bytes).map_err(|e| Error::Sbml(format!("not UTF-8: {e}")))?
    };
    read_str(&text)
}

/// Reads an SBML document from text.
pub fn read_str(xml: &str) -> Result<Document> {
    let (root, packages) = dom::parse(xml)?;
    if root.name != "sbml" {
        return Err(Error::Sbml(format!("root element is <{}>, expected <sbml>", root.name)));
    }
    let level: u32 = parse_num(root.attr("level"), "level")?;
    let version: u32 = parse_num(root.attr("version"), "version")?;
    let model = root.child("model").ok_or_else(|| Error::Sbml("no <model> element".into()))?;
    let model = read_model(model, level)?;
    Ok(Document { level, version, packages, model })
}

fn parse_num<T: std::str::FromStr>(v: Option<&str>, what: &str) -> Result<T> {
    v.ok_or_else(|| Error::Sbml(format!("missing attribute {what}")))?
        .parse()
        .map_err(|_| Error::Sbml(format!("bad value for {what}")))
}

fn opt_f64(node: &Node, key: &str) -> Result<Option<f64>> {
    match node.attr(key) {
        None => Ok(None),
        Some(text) => parse_f64(text).map(Some).map_err(|_| {
            Error::Sbml(format!("bad number `{text}` in attribute {key} of <{}>", node.name))
        }),
    }
}

fn parse_f64(text: &str) -> std::result::Result<f64, ()> {
    match text.trim() {
        "INF" | "inf" | "Infinity" | "+INF" => Ok(f64::INFINITY),
        "-INF" | "-inf" | "-Infinity" => Ok(f64::NEG_INFINITY),
        "NaN" | "nan" => Ok(f64::NAN),
        t => t.parse().map_err(|_| ()),
    }
}

fn flag(node: &Node, key: &str, default: bool) -> bool {
    match node.attr(key) {
        Some("true") | Some("1") => true,
        Some("false") | Some("0") => false,
        _ => default,
    }
}

fn text_attr(node: &Node, key: &str) -> String {
    node.attr(key).unwrap_or_default().to_string()
}

fn list<'a>(parent: &'a Node, list: &str, item: &'a str) -> Vec<&'a Node> {
    parent.child(list).map(|l| l.children_named(item).collect()).unwrap_or_default()
}

fn math_of(node: &Node) -> Result<Option<Expr>> {
    match node.child("math") {
        Some(math) if !math.children.is_empty() => Ok(Some(read_math(math)?)),
        _ => Ok(None),
    }
}

fn annotations(node: &Node) -> Vec<String> {
    let mut out = Vec::new();
    if let Some(annotation) = node.child("annotation") {
        collect_resources(annotation, &mut out);
    }
    out
}

fn collect_resources(node: &Node, out: &mut Vec<String>) {
    if node.name == "rdf:li"
        && let Some(res) = node.attr("rdf:resource")
    {
        out.push(res.to_string());
    }
    for child in &node.children {
        collect_resources(child, out);
    }
}

fn read_model(node: &Node, level: u32) -> Result<SbmlModel> {
    let mut model = SbmlModel {
        id: text_attr(node, "id"),
        name: text_attr(node, "name"),
        conversion_factor: node.attr("conversionFactor").map(str::to_string),
        ..SbmlModel::default()
    };
    for f in list(node, "listOfFunctionDefinitions", "functionDefinition") {
        // Level 3 Version 2 allows a function without math; nothing may call it.
        let Some(math) = f.child("math") else { continue };
        let (args, body) = read_lambda(math)?;
        model.functions.push(FunctionDefinition { id: text_attr(f, "id"), args, body });
    }
    for c in list(node, "listOfCompartments", "compartment") {
        // Level 2 defaults: size 1 when absent is applied by the kinetic compiler.
        model.compartments.push(Compartment {
            id: text_attr(c, "id"),
            name: text_attr(c, "name"),
            size: match opt_f64(c, "size")? {
                Some(v) => Some(v),
                None => opt_f64(c, "volume")?,
            },
            constant: flag(c, "constant", level < 3),
            spatial_dimensions: opt_f64(c, "spatialDimensions")?,
        });
    }
    for s in list(node, "listOfSpecies", "species") {
        model.species.push(Species {
            id: text_attr(s, "id"),
            name: text_attr(s, "name"),
            compartment: text_attr(s, "compartment"),
            initial_amount: opt_f64(s, "initialAmount")?,
            initial_concentration: opt_f64(s, "initialConcentration")?,
            has_only_substance_units: flag(s, "hasOnlySubstanceUnits", false),
            boundary_condition: flag(s, "boundaryCondition", false),
            constant: flag(s, "constant", false),
            conversion_factor: s.attr("conversionFactor").map(str::to_string),
            formula: s.attr("fbc:chemicalFormula").map(str::to_string),
            charge: s.attr("fbc:charge").and_then(|c| c.parse().ok()),
            annotations: annotations(s),
        });
    }
    for p in list(node, "listOfParameters", "parameter") {
        model.parameters.push(read_parameter(p, level)?);
    }
    for ia in list(node, "listOfInitialAssignments", "initialAssignment") {
        if let Some(math) = math_of(ia)? {
            model.initial_assignments.push(InitialAssignment { symbol: text_attr(ia, "symbol"), math });
        }
    }
    if let Some(rules) = node.child("listOfRules") {
        for r in &rules.children {
            let Some(math) = math_of(r)? else { continue };
            let variable = text_attr(r, "variable");
            match r.name.as_str() {
                "assignmentRule" => model.rules.push(Rule::Assignment { variable, math }),
                "rateRule" => model.rules.push(Rule::Rate { variable, math }),
                "algebraicRule" => model.rules.push(Rule::Algebraic { math }),
                other => return Err(Error::Sbml(format!("unsupported rule <{other}>"))),
            }
        }
    }
    for c in list(node, "listOfConstraints", "constraint") {
        if let Some(math) = math_of(c)? {
            model.constraints.push(math);
        }
    }
    for r in list(node, "listOfReactions", "reaction") {
        model.reactions.push(read_reaction(r, level)?);
    }
    for e in list(node, "listOfEvents", "event") {
        model.events.push(read_event(e, level)?);
    }
    read_fbc(node, &mut model)?;
    read_groups(node, &mut model);
    Ok(model)
}

fn read_parameter(p: &Node, level: u32) -> Result<Parameter> {
    Ok(Parameter {
        id: text_attr(p, "id"),
        name: text_attr(p, "name"),
        value: opt_f64(p, "value")?,
        constant: flag(p, "constant", level < 3 || p.name == "localParameter"),
    })
}

fn read_species_ref(s: &Node) -> Result<SpeciesReference> {
    Ok(SpeciesReference {
        species: text_attr(s, "species"),
        id: s.attr("id").map(str::to_string),
        stoichiometry: opt_f64(s, "stoichiometry")?,
        constant: flag(s, "constant", true),
        stoichiometry_math: match s.child("stoichiometryMath") {
            Some(m) => math_of(m)?,
            None => None,
        },
    })
}

fn read_reaction(r: &Node, level: u32) -> Result<Reaction> {
    let mut out = Reaction {
        id: text_attr(r, "id"),
        name: text_attr(r, "name"),
        reversible: flag(r, "reversible", true),
        fast: flag(r, "fast", false),
        lower_bound: r.attr("fbc:lowerFluxBound").map(str::to_string),
        upper_bound: r.attr("fbc:upperFluxBound").map(str::to_string),
        annotations: annotations(r),
        ..Reaction::default()
    };
    for s in list(r, "listOfReactants", "speciesReference") {
        out.reactants.push(read_species_ref(s)?);
    }
    for s in list(r, "listOfProducts", "speciesReference") {
        out.products.push(read_species_ref(s)?);
    }
    for s in list(r, "listOfModifiers", "modifierSpeciesReference") {
        out.modifiers.push(text_attr(s, "species"));
    }
    if let Some(law) = r.child("kineticLaw") {
        out.kinetic_law = math_of(law)?;
        for p in list(law, "listOfLocalParameters", "localParameter") {
            out.local_parameters.push(read_parameter(p, level)?);
        }
        for p in list(law, "listOfParameters", "parameter") {
            out.local_parameters.push(read_parameter(p, level)?);
        }
    }
    if let Some(gpa) = r.child("fbc:geneProductAssociation")
        && let Some(inner) = gpa.children.first()
    {
        out.gene_rule = Some(gpa_rule(inner));
    }
    if out.gene_rule.is_none()
        && let Some(notes) = r.child("notes")
    {
        out.gene_rule = legacy_gene_rule(&notes.deep_text());
    }
    Ok(out)
}

/// Text rule from nested `fbc:and`, `fbc:or` and `fbc:geneProductRef`.
fn gpa_rule(node: &Node) -> String {
    match node.name.as_str() {
        "fbc:geneProductRef" => text_attr(node, "fbc:geneProduct"),
        "fbc:and" | "fbc:or" => {
            let sep = if node.name == "fbc:and" { " and " } else { " or " };
            let parts: Vec<String> = node.children.iter().map(|c| format!("({})", gpa_rule(c))).collect();
            parts.join(sep)
        }
        _ => String::new(),
    }
}

/// COBRA-style `GENE_ASSOCIATION: ...` or `GPR: ...` line in the reaction notes.
fn legacy_gene_rule(notes: &str) -> Option<String> {
    for key in ["GENE_ASSOCIATION:", "GENE ASSOCIATION:", "GPR:"] {
        if let Some(pos) = notes.find(key) {
            let rest = &notes[pos + key.len()..];
            let line = rest.lines().next().unwrap_or_default().trim();
            if !line.is_empty() {
                return Some(line.to_string());
            }
        }
    }
    None
}

fn read_event(e: &Node, level: u32) -> Result<Event> {
    // Level 3 Version 2 allows an event without a trigger or a trigger without math; neither fires.
    let trigger_node = e.child("trigger");
    let trigger = match trigger_node {
        Some(t) => math_of(t)?.unwrap_or(Expr::Num(0.0)),
        None => Expr::Num(0.0),
    };
    let delay = match e.child("delay") {
        Some(d) => math_of(d)?,
        None => None,
    };
    let priority = match e.child("priority") {
        Some(p) => math_of(p)?,
        None => None,
    };
    let mut assignments = Vec::new();
    for ea in list(e, "listOfEventAssignments", "eventAssignment") {
        if let Some(math) = math_of(ea)? {
            assignments.push(EventAssignment { variable: text_attr(ea, "variable"), math });
        }
    }
    Ok(Event {
        id: text_attr(e, "id"),
        trigger,
        initial_value: trigger_node.is_none_or(|t| flag(t, "initialValue", level < 3)),
        persistent: trigger_node.is_none_or(|t| flag(t, "persistent", level < 3)),
        delay,
        priority,
        use_values_from_trigger_time: flag(e, "useValuesFromTriggerTime", true),
        assignments,
    })
}

fn read_fbc(node: &Node, model: &mut SbmlModel) -> Result<()> {
    if let Some(objs) = node.child("fbc:listOfObjectives") {
        model.active_objective = objs.attr("fbc:activeObjective").map(str::to_string);
        for o in objs.children_named("fbc:objective") {
            let mut objective = Objective {
                id: text_attr(o, "fbc:id"),
                maximize: o.attr("fbc:type").is_some_and(|t| t.eq_ignore_ascii_case("maximize")),
                flux_objectives: Vec::new(),
            };
            if let Some(fluxes) = o.child("fbc:listOfFluxObjectives") {
                for f in fluxes.children_named("fbc:fluxObjective") {
                    let coef = match opt_f64(f, "fbc:coefficient")? {
                        Some(c) => c,
                        None => 1.0,
                    };
                    objective.flux_objectives.push((text_attr(f, "fbc:reaction"), coef));
                }
            }
            model.objectives.push(objective);
        }
    }
    if let Some(genes) = node.child("fbc:listOfGeneProducts") {
        for g in genes.children_named("fbc:geneProduct") {
            model.gene_products.push(GeneProduct {
                id: text_attr(g, "fbc:id"),
                label: text_attr(g, "fbc:label"),
                name: text_attr(g, "fbc:name"),
            });
        }
    }
    // FBC version 1 stores bounds as a separate list of flux bounds with literal values.
    if let Some(bounds) = node.child("fbc:listOfFluxBounds") {
        for b in bounds.children_named("fbc:fluxBound") {
            let rid = text_attr(b, "fbc:reaction");
            let op = text_attr(b, "fbc:operation");
            let value = opt_f64(b, "fbc:value")?.unwrap_or(f64::NAN);
            let pid = format!("__fbc_bound_{rid}_{op}");
            model.parameters.push(Parameter {
                id: pid.clone(),
                name: String::new(),
                value: Some(value),
                constant: true,
            });
            if let Some(r) = model.reactions.iter_mut().find(|r| r.id == rid) {
                match op.as_str() {
                    "lessEqual" => r.upper_bound = Some(pid),
                    "greaterEqual" => r.lower_bound = Some(pid),
                    "equal" => {
                        r.lower_bound = Some(pid.clone());
                        r.upper_bound = Some(pid);
                    }
                    _ => {}
                }
            }
        }
    }
    Ok(())
}

fn read_groups(node: &Node, model: &mut SbmlModel) {
    if let Some(groups) = node.child("groups:listOfGroups") {
        for g in groups.children_named("groups:group") {
            let mut group = Group {
                id: text_attr(g, "groups:id"),
                name: text_attr(g, "groups:name"),
                kind: text_attr(g, "groups:kind"),
                members: Vec::new(),
            };
            if let Some(members) = g.child("groups:listOfMembers") {
                for m in members.children_named("groups:member") {
                    if let Some(id) = m.attr("groups:idRef") {
                        group.members.push(id.to_string());
                    }
                }
            }
            model.groups.push(group);
        }
    }
}
