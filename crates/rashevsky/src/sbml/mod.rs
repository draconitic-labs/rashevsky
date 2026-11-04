//! SBML reading and writing: core Levels 2 and 3, with the FBC (versions 1 to 3) and groups packages.

pub mod dom;
pub mod mathml;
pub mod read;
pub mod write;

use crate::kinetic::expr::Expr;
use serde::{Deserialize, Serialize};

/// A parsed SBML document.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Document {
    /// SBML level.
    pub level: u32,
    /// SBML version.
    pub version: u32,
    /// Packages declared on the root element.
    pub packages: Vec<String>,
    /// The model.
    pub model: SbmlModel,
}

/// The model element of an SBML document.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SbmlModel {
    /// Model identifier.
    pub id: String,
    /// Human-readable name.
    pub name: String,
    /// Model-wide conversion factor parameter.
    pub conversion_factor: Option<String>,
    /// Function definitions.
    pub functions: Vec<FunctionDefinition>,
    /// Compartments.
    pub compartments: Vec<Compartment>,
    /// Species.
    pub species: Vec<Species>,
    /// Global parameters.
    pub parameters: Vec<Parameter>,
    /// Initial assignments.
    pub initial_assignments: Vec<InitialAssignment>,
    /// Rules, in document order.
    pub rules: Vec<Rule>,
    /// Constraints; they are read, never enforced during simulation.
    pub constraints: Vec<Expr>,
    /// Reactions.
    pub reactions: Vec<Reaction>,
    /// Events.
    pub events: Vec<Event>,
    /// FBC objectives.
    pub objectives: Vec<Objective>,
    /// Identifier of the active FBC objective.
    pub active_objective: Option<String>,
    /// FBC gene products.
    pub gene_products: Vec<GeneProduct>,
    /// Groups-package groups.
    pub groups: Vec<Group>,
}

/// A function definition.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FunctionDefinition {
    /// Identifier.
    pub id: String,
    /// Bound variable names.
    pub args: Vec<String>,
    /// Lambda body.
    pub body: Expr,
}

/// A compartment.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Compartment {
    /// Identifier.
    pub id: String,
    /// Name.
    pub name: String,
    /// Size, when given.
    pub size: Option<f64>,
    /// Whether the size is constant.
    pub constant: bool,
    /// Spatial dimensions, when given.
    pub spatial_dimensions: Option<f64>,
}

/// A species.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Species {
    /// Identifier.
    pub id: String,
    /// Name.
    pub name: String,
    /// Compartment identifier.
    pub compartment: String,
    /// Initial amount, when given.
    pub initial_amount: Option<f64>,
    /// Initial concentration, when given.
    pub initial_concentration: Option<f64>,
    /// Whether math refers to the amount rather than the concentration.
    pub has_only_substance_units: bool,
    /// Whether the species is outside the reaction system's balance.
    pub boundary_condition: bool,
    /// Whether the value never changes.
    pub constant: bool,
    /// Species conversion factor.
    pub conversion_factor: Option<String>,
    /// FBC formula.
    pub formula: Option<String>,
    /// FBC charge.
    pub charge: Option<i32>,
    /// `bqbiol:is` annotation resources.
    pub annotations: Vec<String>,
}

/// A parameter.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Parameter {
    /// Identifier.
    pub id: String,
    /// Name.
    pub name: String,
    /// Value, when given.
    pub value: Option<f64>,
    /// Whether the value is constant.
    pub constant: bool,
}

/// An initial assignment.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct InitialAssignment {
    /// Symbol assigned.
    pub symbol: String,
    /// Value.
    pub math: Expr,
}

/// A rule.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Rule {
    /// `variable = math`.
    Assignment {
        /// Assigned symbol.
        variable: String,
        /// Expression.
        math: Expr,
    },
    /// `d(variable)/dt = math`.
    Rate {
        /// Symbol.
        variable: String,
        /// Expression.
        math: Expr,
    },
    /// `0 = math`.
    Algebraic {
        /// Expression.
        math: Expr,
    },
}

/// A species reference in a reaction.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SpeciesReference {
    /// Referenced species.
    pub species: String,
    /// Identifier of the reference, usable as a symbol in math.
    pub id: Option<String>,
    /// Stoichiometry, when given.
    pub stoichiometry: Option<f64>,
    /// Whether the stoichiometry is constant.
    pub constant: bool,
    /// Level 2 stoichiometry expression.
    pub stoichiometry_math: Option<Expr>,
}

/// A reaction.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Reaction {
    /// Identifier.
    pub id: String,
    /// Name.
    pub name: String,
    /// Reversibility flag.
    pub reversible: bool,
    /// Fast flag (removed in Level 3 Version 2).
    pub fast: bool,
    /// Reactants.
    pub reactants: Vec<SpeciesReference>,
    /// Products.
    pub products: Vec<SpeciesReference>,
    /// Modifier species.
    pub modifiers: Vec<String>,
    /// Kinetic law expression.
    pub kinetic_law: Option<Expr>,
    /// Local parameters of the kinetic law.
    pub local_parameters: Vec<Parameter>,
    /// FBC lower bound parameter.
    pub lower_bound: Option<String>,
    /// FBC upper bound parameter.
    pub upper_bound: Option<String>,
    /// Gene association as text rule over gene product identifiers.
    pub gene_rule: Option<String>,
    /// `bqbiol:is` annotation resources.
    pub annotations: Vec<String>,
}

/// One assignment of an event.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EventAssignment {
    /// Symbol assigned.
    pub variable: String,
    /// Value.
    pub math: Expr,
}

/// An event.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Event {
    /// Identifier.
    pub id: String,
    /// Trigger condition.
    pub trigger: Expr,
    /// Whether the trigger is true before time zero.
    pub initial_value: bool,
    /// Whether a pending execution cancels when the trigger goes false.
    pub persistent: bool,
    /// Delay expression.
    pub delay: Option<Expr>,
    /// Priority expression.
    pub priority: Option<Expr>,
    /// Whether assignments are computed at trigger time.
    pub use_values_from_trigger_time: bool,
    /// Assignments.
    pub assignments: Vec<EventAssignment>,
}

/// An FBC objective.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Objective {
    /// Identifier.
    pub id: String,
    /// True for maximize.
    pub maximize: bool,
    /// `(reaction, coefficient)` pairs.
    pub flux_objectives: Vec<(String, f64)>,
}

/// An FBC gene product.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct GeneProduct {
    /// SBML identifier.
    pub id: String,
    /// `fbc:label`, the gene's external identifier.
    pub label: String,
    /// Name.
    pub name: String,
}

/// A group from the groups package.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Group {
    /// Identifier.
    pub id: String,
    /// Name.
    pub name: String,
    /// `partonomy`, `classification` or `collection`.
    pub kind: String,
    /// Identifiers of the members.
    pub members: Vec<String>,
}
