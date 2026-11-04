//! SBML writer: Level 3 Version 1 core with FBC version 2 and groups version 1.

use super::Document;
use super::mathml::{write_condition, write_expr};
use crate::kinetic::expr::Expr;
use std::fmt::Write as _;

/// Escapes text for use in XML attribute values and content.
pub fn escape(text: &str) -> String {
    text.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

fn num(v: f64) -> String {
    if v.is_infinite() {
        if v > 0.0 { "INF".into() } else { "-INF".into() }
    } else if v.is_nan() {
        "NaN".into()
    } else {
        format!("{v:?}")
    }
}

const MATH: &str = "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">";

fn math(e: &Expr, out: &mut String) {
    out.push_str(MATH);
    write_expr(e, out);
    out.push_str("</math>");
}

fn condition(e: &Expr, out: &mut String) {
    out.push_str(MATH);
    write_condition(e, out);
    out.push_str("</math>");
}

/// `bqbiol:is` resources as an annotation on the element with this `metaid`.
fn annotation(metaid: &str, resources: &[String], out: &mut String) {
    let _ = write!(
        out,
        "<annotation><rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\"><rdf:Description rdf:about=\"#{metaid}\"><bqbiol:is><rdf:Bag>"
    );
    for r in resources {
        let _ = write!(out, "<rdf:li rdf:resource=\"{}\"/>", escape(r));
    }
    out.push_str("</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>");
}

fn lambda(args: &[String], body: &Expr, out: &mut String) {
    out.push_str(MATH);
    out.push_str("<lambda>");
    for arg in args {
        let _ = write!(out, "<bvar><ci>{}</ci></bvar>", escape(arg));
    }
    write_expr(body, out);
    out.push_str("</lambda></math>");
}

/// Renders a document as SBML Level 3 Version 1 text, keeping every object the reader understands.
pub fn write_string(doc: &Document) -> String {
    let m = &doc.model;
    let older_level = doc.level < 3;
    // A Level 2 stoichiometry expression becomes a rule on a named species reference.
    let mut stoichiometry_rules: Vec<(String, &Expr)> = Vec::new();
    let mut reference_ids: std::collections::HashMap<(usize, bool, usize), String> = std::collections::HashMap::new();
    for (j, r) in m.reactions.iter().enumerate() {
        for (is_product, refs) in [(false, &r.reactants), (true, &r.products)] {
            for (k, sr) in refs.iter().enumerate() {
                if let Some(e) = &sr.stoichiometry_math {
                    let id = sr.id.clone().unwrap_or_else(|| format!("stoichiometry_{j}_{}_{k}", u8::from(is_product)));
                    reference_ids.insert((j, is_product, k), id.clone());
                    stoichiometry_rules.push((id, e));
                }
            }
        }
    }
    let mut s = String::with_capacity(1 << 16);
    s.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    // Only the packages a model uses are declared; some readers stop at an unused package.
    let uses_fbc = !m.objectives.is_empty()
        || !m.gene_products.is_empty()
        || m.reactions.iter().any(|r| r.lower_bound.is_some() || r.upper_bound.is_some() || r.gene_rule.is_some())
        || m.species.iter().any(|sp| sp.formula.is_some() || sp.charge.is_some());
    let uses_groups = !m.groups.is_empty();
    s.push_str("<sbml xmlns=\"http://www.sbml.org/sbml/level3/version1/core\" level=\"3\" version=\"1\"");
    if uses_fbc {
        s.push_str(" xmlns:fbc=\"http://www.sbml.org/sbml/level3/version1/fbc/version2\" fbc:required=\"false\"");
    }
    if uses_groups {
        s.push_str(" xmlns:groups=\"http://www.sbml.org/sbml/level3/version1/groups/version1\" groups:required=\"false\"");
    }
    s.push_str(">\n");
    let _ = write!(s, "<model id=\"{}\"", escape(&m.id));
    if uses_fbc {
        s.push_str(" fbc:strict=\"false\"");
    }
    if !m.name.is_empty() {
        let _ = write!(s, " name=\"{}\"", escape(&m.name));
    }
    if let Some(cf) = &m.conversion_factor {
        let _ = write!(s, " conversionFactor=\"{}\"", escape(cf));
    }
    s.push_str(">\n");
    if !m.functions.is_empty() {
        s.push_str("<listOfFunctionDefinitions>\n");
        for f in &m.functions {
            let _ = write!(s, "<functionDefinition id=\"{}\">", escape(&f.id));
            lambda(&f.args, &f.body, &mut s);
            s.push_str("</functionDefinition>\n");
        }
        s.push_str("</listOfFunctionDefinitions>\n");
    }
    if !m.compartments.is_empty() {
        s.push_str("<listOfCompartments>\n");
        for c in &m.compartments {
            let _ = write!(s, "<compartment id=\"{}\" constant=\"{}\"", escape(&c.id), c.constant);
            if !c.name.is_empty() {
                let _ = write!(s, " name=\"{}\"", escape(&c.name));
            }
            // Level 2 gives a missing size the default of one, and Level 3 leaves it undefined.
            if let Some(v) = c.size.or(older_level.then_some(1.0)) {
                let _ = write!(s, " size=\"{}\"", num(v));
            }
            if let Some(v) = c.spatial_dimensions.or(older_level.then_some(3.0)) {
                let _ = write!(s, " spatialDimensions=\"{}\"", num(v));
            }
            s.push_str("/>\n");
        }
        s.push_str("</listOfCompartments>\n");
    }
    if !m.species.is_empty() {
        s.push_str("<listOfSpecies>\n");
        for (i, sp) in m.species.iter().enumerate() {
            let _ = write!(
                s,
                "<species id=\"{}\" compartment=\"{}\" hasOnlySubstanceUnits=\"{}\" boundaryCondition=\"{}\" constant=\"{}\"",
                escape(&sp.id),
                escape(&sp.compartment),
                sp.has_only_substance_units,
                sp.boundary_condition,
                sp.constant
            );
            if !sp.name.is_empty() {
                let _ = write!(s, " name=\"{}\"", escape(&sp.name));
            }
            if let Some(v) = sp.initial_amount {
                let _ = write!(s, " initialAmount=\"{}\"", num(v));
            }
            if let Some(v) = sp.initial_concentration {
                let _ = write!(s, " initialConcentration=\"{}\"", num(v));
            }
            if let Some(cf) = &sp.conversion_factor {
                let _ = write!(s, " conversionFactor=\"{}\"", escape(cf));
            }
            if let Some(f) = &sp.formula {
                let _ = write!(s, " fbc:chemicalFormula=\"{}\"", escape(f));
            }
            if let Some(c) = sp.charge {
                let _ = write!(s, " fbc:charge=\"{c}\"");
            }
            if sp.annotations.is_empty() {
                s.push_str("/>\n");
            } else {
                let metaid = format!("meta_species_{i}");
                let _ = write!(s, " metaid=\"{metaid}\">");
                annotation(&metaid, &sp.annotations, &mut s);
                s.push_str("</species>\n");
            }
        }
        s.push_str("</listOfSpecies>\n");
    }
    if !m.parameters.is_empty() {
        s.push_str("<listOfParameters>\n");
        for p in &m.parameters {
            let _ = write!(s, "<parameter id=\"{}\" constant=\"{}\"", escape(&p.id), p.constant);
            if let Some(v) = p.value {
                let _ = write!(s, " value=\"{}\"", num(v));
            }
            s.push_str("/>\n");
        }
        s.push_str("</listOfParameters>\n");
    }
    if !m.initial_assignments.is_empty() {
        s.push_str("<listOfInitialAssignments>\n");
        for ia in &m.initial_assignments {
            let _ = write!(s, "<initialAssignment symbol=\"{}\">", escape(&ia.symbol));
            math(&ia.math, &mut s);
            s.push_str("</initialAssignment>\n");
        }
        s.push_str("</listOfInitialAssignments>\n");
    }
    if !m.rules.is_empty() || !stoichiometry_rules.is_empty() {
        s.push_str("<listOfRules>\n");
        for (id, e) in &stoichiometry_rules {
            let _ = write!(s, "<assignmentRule variable=\"{}\">", escape(id));
            math(e, &mut s);
            s.push_str("</assignmentRule>\n");
        }
        for rule in &m.rules {
            match rule {
                super::Rule::Assignment { variable, math: e } => {
                    let _ = write!(s, "<assignmentRule variable=\"{}\">", escape(variable));
                    math(e, &mut s);
                    s.push_str("</assignmentRule>\n");
                }
                super::Rule::Rate { variable, math: e } => {
                    let _ = write!(s, "<rateRule variable=\"{}\">", escape(variable));
                    math(e, &mut s);
                    s.push_str("</rateRule>\n");
                }
                super::Rule::Algebraic { math: e } => {
                    s.push_str("<algebraicRule>");
                    math(e, &mut s);
                    s.push_str("</algebraicRule>\n");
                }
            }
        }
        s.push_str("</listOfRules>\n");
    }
    if !m.reactions.is_empty() {
        s.push_str("<listOfReactions>\n");
        for (j, r) in m.reactions.iter().enumerate() {
            let _ = write!(s, "<reaction id=\"{}\" reversible=\"{}\" fast=\"{}\"", escape(&r.id), r.reversible, r.fast);
            if !r.annotations.is_empty() {
                let _ = write!(s, " metaid=\"meta_reaction_{j}\"");
            }
            if !r.name.is_empty() {
                let _ = write!(s, " name=\"{}\"", escape(&r.name));
            }
            if let Some(b) = &r.lower_bound {
                let _ = write!(s, " fbc:lowerFluxBound=\"{}\"", escape(b));
            }
            if let Some(b) = &r.upper_bound {
                let _ = write!(s, " fbc:upperFluxBound=\"{}\"", escape(b));
            }
            s.push_str(">\n");
            if !r.annotations.is_empty() {
                annotation(&format!("meta_reaction_{j}"), &r.annotations, &mut s);
            }
            for (tag, is_product, refs) in [("listOfReactants", false, &r.reactants), ("listOfProducts", true, &r.products)] {
                if refs.is_empty() {
                    continue;
                }
                let _ = writeln!(s, "<{tag}>");
                for (k, sr) in refs.iter().enumerate() {
                    let named = reference_ids.get(&(j, is_product, k));
                    let _ = write!(
                        s,
                        "<speciesReference species=\"{}\" constant=\"{}\"",
                        escape(&sr.species),
                        sr.constant && named.is_none()
                    );
                    if let Some(id) = named.or(sr.id.as_ref()) {
                        let _ = write!(s, " id=\"{}\"", escape(id));
                    }
                    // Level 2 gives a missing stoichiometry the default of one.
                    if named.is_none()
                        && let Some(v) = sr.stoichiometry.or(older_level.then_some(1.0))
                    {
                        let _ = write!(s, " stoichiometry=\"{}\"", num(v));
                    }
                    s.push_str("/>\n");
                }
                let _ = writeln!(s, "</{tag}>");
            }
            if !r.modifiers.is_empty() {
                s.push_str("<listOfModifiers>\n");
                for sp in &r.modifiers {
                    let _ = writeln!(s, "<modifierSpeciesReference species=\"{}\"/>", escape(sp));
                }
                s.push_str("</listOfModifiers>\n");
            }
            if let Some(law) = &r.kinetic_law {
                s.push_str("<kineticLaw>");
                math(law, &mut s);
                if !r.local_parameters.is_empty() {
                    s.push_str("<listOfLocalParameters>");
                    for p in &r.local_parameters {
                        let _ = write!(s, "<localParameter id=\"{}\"", escape(&p.id));
                        if let Some(v) = p.value {
                            let _ = write!(s, " value=\"{}\"", num(v));
                        }
                        s.push_str("/>");
                    }
                    s.push_str("</listOfLocalParameters>");
                }
                s.push_str("</kineticLaw>\n");
            }
            if let Some(rule) = &r.gene_rule
                && !rule.trim().is_empty()
            {
                write_gene_rule(rule, &mut s);
            }
            s.push_str("</reaction>\n");
        }
        s.push_str("</listOfReactions>\n");
    }
    if !m.events.is_empty() {
        s.push_str("<listOfEvents>\n");
        for e in &m.events {
            let _ = write!(
                s,
                "<event id=\"{}\" useValuesFromTriggerTime=\"{}\">",
                escape(&e.id),
                e.use_values_from_trigger_time
            );
            let _ = write!(s, "<trigger initialValue=\"{}\" persistent=\"{}\">", e.initial_value, e.persistent);
            condition(&e.trigger, &mut s);
            s.push_str("</trigger>");
            if let Some(p) = &e.priority {
                s.push_str("<priority>");
                math(p, &mut s);
                s.push_str("</priority>");
            }
            if let Some(d) = &e.delay {
                s.push_str("<delay>");
                math(d, &mut s);
                s.push_str("</delay>");
            }
            s.push_str("<listOfEventAssignments>");
            for a in &e.assignments {
                let _ = write!(s, "<eventAssignment variable=\"{}\">", escape(&a.variable));
                math(&a.math, &mut s);
                s.push_str("</eventAssignment>");
            }
            s.push_str("</listOfEventAssignments></event>\n");
        }
        s.push_str("</listOfEvents>\n");
    }
    if !m.objectives.is_empty() {
        let active = m.active_objective.clone().unwrap_or_else(|| m.objectives[0].id.clone());
        let _ = writeln!(s, "<fbc:listOfObjectives fbc:activeObjective=\"{}\">", escape(&active));
        for o in &m.objectives {
            let _ = writeln!(
                s,
                "<fbc:objective fbc:id=\"{}\" fbc:type=\"{}\">",
                escape(&o.id),
                if o.maximize { "maximize" } else { "minimize" }
            );
            s.push_str("<fbc:listOfFluxObjectives>\n");
            for (rid, coef) in &o.flux_objectives {
                let _ = writeln!(
                    s,
                    "<fbc:fluxObjective fbc:reaction=\"{}\" fbc:coefficient=\"{}\"/>",
                    escape(rid),
                    num(*coef)
                );
            }
            s.push_str("</fbc:listOfFluxObjectives>\n</fbc:objective>\n");
        }
        s.push_str("</fbc:listOfObjectives>\n");
    }
    if !m.gene_products.is_empty() {
        s.push_str("<fbc:listOfGeneProducts>\n");
        for g in &m.gene_products {
            let _ = write!(s, "<fbc:geneProduct fbc:id=\"{}\" fbc:label=\"{}\"", escape(&g.id), escape(&g.label));
            if !g.name.is_empty() {
                let _ = write!(s, " fbc:name=\"{}\"", escape(&g.name));
            }
            s.push_str("/>\n");
        }
        s.push_str("</fbc:listOfGeneProducts>\n");
    }
    if !m.groups.is_empty() {
        s.push_str("<groups:listOfGroups>\n");
        for g in &m.groups {
            let _ = write!(s, "<groups:group groups:id=\"{}\" groups:kind=\"{}\"", escape(&g.id), escape(&g.kind));
            if !g.name.is_empty() {
                let _ = write!(s, " groups:name=\"{}\"", escape(&g.name));
            }
            s.push_str("><groups:listOfMembers>\n");
            for member in &g.members {
                let _ = writeln!(s, "<groups:member groups:idRef=\"{}\"/>", escape(member));
            }
            s.push_str("</groups:listOfMembers></groups:group>\n");
        }
        s.push_str("</groups:listOfGroups>\n");
    }
    s.push_str("</model>\n</sbml>\n");
    s
}

fn write_gene_rule(rule: &str, s: &mut String) {
    let mut idx = crate::gpr::GeneIndex::default();
    let Ok(gpr) = crate::gpr::Gpr::parse(rule, &mut |g| idx.intern(g)) else {
        return;
    };
    s.push_str("<fbc:geneProductAssociation>");
    write_gpr(&gpr, &idx, s);
    s.push_str("</fbc:geneProductAssociation>\n");
}

fn write_gpr(gpr: &crate::gpr::Gpr, idx: &crate::gpr::GeneIndex, s: &mut String) {
    use crate::gpr::Gpr;
    match gpr {
        Gpr::None => {}
        Gpr::Gene(g) => {
            let _ = write!(s, "<fbc:geneProductRef fbc:geneProduct=\"{}\"/>", escape(&idx.ids[*g]));
        }
        Gpr::And(items) => {
            s.push_str("<fbc:and>");
            items.iter().for_each(|i| write_gpr(i, idx, s));
            s.push_str("</fbc:and>");
        }
        Gpr::Or(items) => {
            s.push_str("<fbc:or>");
            items.iter().for_each(|i| write_gpr(i, idx, s));
            s.push_str("</fbc:or>");
        }
    }
}
