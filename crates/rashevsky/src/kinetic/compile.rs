//! Compiles an SBML model into a simulatable form.
//!
//! Every symbol gets a slot in a value array. Each species has two slots: `S` is the value math sees
//! (a concentration, or an amount with `hasOnlySubstanceUnits`) and `S#a` its amount. The integrated
//! state consists of the amounts of reaction-driven species and of every symbol with a rate rule,
//! algebraic rule or event assignment. All other values derive from the state through expressions
//! evaluated in dependency order (assignment rules, species views, reaction rates, `rateOf`).

use super::diff::{self, Var, add, div, mul, neg, sub};
use super::expr::{Expr, Op, apply_op};
use crate::error::{Error, Result};
use crate::sbml::{self, Document, Rule, SbmlModel};
use std::collections::{HashMap, HashSet};

/// Avogadro's constant as fixed by SBML Level 3.
pub const AVOGADRO: f64 = 6.02214179e23;

/// An expression with symbols resolved to slots.
#[derive(Clone, Debug)]
pub enum C {
    /// Constant.
    Num(f64),
    /// Value of a slot.
    Slot(usize),
    /// Simulation time.
    Time,
    /// Built-in operator.
    Ap(Op, Vec<C>),
    /// Piecewise function, evaluated lazily.
    Piece(Vec<(C, C)>, Option<Box<C>>),
    /// Value of an expression a lag ago.
    Delay(Box<C>, Box<C>),
}

/// Past values for `delay`.
pub trait DelayContext {
    /// Evaluates `expr` at time `t`, which may lie in the past.
    fn value_at(&self, expr: &C, t: f64) -> f64;
}

/// A context for models without delays.
pub struct NoDelay;

impl DelayContext for NoDelay {
    fn value_at(&self, _expr: &C, _t: f64) -> f64 {
        f64::NAN
    }
}

/// Evaluates a compiled expression.
pub fn eval(c: &C, vals: &[f64], t: f64, ctx: &dyn DelayContext) -> f64 {
    match c {
        C::Num(v) => *v,
        C::Slot(s) => vals[*s],
        C::Time => t,
        C::Ap(op, args) => match (op, args.as_slice()) {
            (Op::Add, [a, b]) => eval(a, vals, t, ctx) + eval(b, vals, t, ctx),
            (Op::Sub, [a, b]) => eval(a, vals, t, ctx) - eval(b, vals, t, ctx),
            (Op::Mul, [a, b]) => eval(a, vals, t, ctx) * eval(b, vals, t, ctx),
            (Op::Div, [a, b]) => eval(a, vals, t, ctx) / eval(b, vals, t, ctx),
            (Op::Neg, [a]) => -eval(a, vals, t, ctx),
            (Op::Add, many) => many.iter().map(|a| eval(a, vals, t, ctx)).sum(),
            (Op::Mul, many) => many.iter().map(|a| eval(a, vals, t, ctx)).product(),
            (op, many) => {
                let mut small = [0.0; 6];
                if many.len() <= 6 {
                    for (slot, a) in small.iter_mut().zip(many) {
                        *slot = eval(a, vals, t, ctx);
                    }
                    apply_op(*op, &small[..many.len()])
                } else {
                    let args: Vec<f64> = many.iter().map(|a| eval(a, vals, t, ctx)).collect();
                    apply_op(*op, &args)
                }
            }
        },
        C::Piece(pieces, otherwise) => {
            for (value, cond) in pieces {
                if eval(cond, vals, t, ctx) != 0.0 {
                    return eval(value, vals, t, ctx);
                }
            }
            otherwise.as_ref().map_or(f64::NAN, |o| eval(o, vals, t, ctx))
        }
        C::Delay(x, lag) => {
            let lag = eval(lag, vals, t, ctx);
            ctx.value_at(x, t - lag)
        }
    }
}

/// How a species is represented.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpeciesKind {
    /// Its amount is a state variable.
    Amount,
    /// Its math value has a rate rule, and that value is the state variable.
    RateRule,
    /// Its math value is determined by an assignment rule.
    Assigned,
    /// Its math value is determined by an algebraic rule.
    Algebraic,
    /// Its math value never changes.
    Constant,
}

/// Slots and representation of a species.
#[derive(Clone, Debug)]
pub struct SpeciesInfo {
    /// Identifier.
    pub id: String,
    /// Slot of the value math sees.
    pub view: usize,
    /// Slot of the amount.
    pub amount: usize,
    /// Slot of the compartment size.
    pub compartment: usize,
    /// Whether math sees the amount.
    pub has_only_substance_units: bool,
    /// Representation.
    pub kind: SpeciesKind,
}

/// A stoichiometric coupling of a reaction to a state variable.
#[derive(Clone, Debug)]
pub struct StoichEntry {
    /// State index.
    pub state: usize,
    /// Reaction index.
    pub reaction: usize,
    /// +1 for products, -1 for reactants.
    pub sign: f64,
    /// Constant coefficient, or the slot of a variable one.
    pub coef: Coef,
    /// Slot of the conversion factor, if any.
    pub conversion: Option<usize>,
}

/// A stoichiometric coefficient.
#[derive(Clone, Copy, Debug)]
pub enum Coef {
    /// Fixed value.
    Const(f64),
    /// Value held in a slot.
    Slot(usize),
}

/// A reaction's identifiers and rate slot.
#[derive(Clone, Debug)]
pub struct ReactionInfo {
    /// Identifier.
    pub id: String,
    /// Slot of the rate.
    pub rate: usize,
    /// Fully inlined rate expression over slot names.
    pub law: Expr,
}

/// Assignment target of an event.
#[derive(Clone, Copy, Debug)]
pub enum Target {
    /// The value is written to the state variable as is.
    State(usize),
    /// The value is in species units; the amount state is the value times the size unless `hosu`.
    Amount {
        /// State index of the amount.
        state: usize,
        /// Slot of the compartment size.
        compartment: usize,
        /// `hasOnlySubstanceUnits` of the species.
        hosu: bool,
    },
}

/// A compiled event.
#[derive(Clone, Debug)]
pub struct CompiledEvent {
    /// Identifier.
    pub id: String,
    /// Trigger.
    pub trigger: C,
    /// Trigger value before the start of the simulation.
    pub initial_value: bool,
    /// Whether execution survives the trigger going false.
    pub persistent: bool,
    /// Delay.
    pub delay: Option<C>,
    /// Priority.
    pub priority: Option<C>,
    /// Compute assignments when triggered rather than when executed.
    pub use_values_from_trigger_time: bool,
    /// Targets with their expressions.
    pub assignments: Vec<(Target, C)>,
}

/// A compiled kinetic model.
#[derive(Clone, Debug)]
pub struct KineticModel {
    /// Model identifier.
    pub id: String,
    /// Model name.
    pub name: String,
    /// Names of all slots.
    pub symbols: Vec<String>,
    /// Slot of a name.
    pub index: HashMap<String, usize>,
    /// Slot of each state variable.
    pub state_slots: Vec<usize>,
    /// Mass-matrix diagonal.
    pub mass: Vec<f64>,
    /// Species.
    pub species: Vec<SpeciesInfo>,
    /// Reactions.
    pub reactions: Vec<ReactionInfo>,
    /// Stoichiometric couplings.
    pub stoich: Vec<StoichEntry>,
    /// `(state, expression)` of each rate rule.
    pub rate_rules: Vec<(usize, C)>,
    /// `(state, residual)` of each algebraic rule.
    pub algebraic: Vec<(usize, C)>,
    /// Events.
    pub events: Vec<CompiledEvent>,
    /// Derived slots in evaluation order.
    pub derived: Vec<(usize, C)>,
    /// All slots at the initial time.
    pub initial: Vec<f64>,
    /// Right-hand side per state as an expression over state and constant slot names.
    pub flat_rhs: Vec<Expr>,
    /// Sparse Jacobian entries `(row, column, expression)`, when analytic differentiation succeeded.
    pub jacobian: Option<Vec<(usize, usize, C)>>,
    /// Time derivative entries.
    pub time_derivative: Vec<(usize, C)>,
    /// Whether any expression uses `delay`.
    pub has_delay: bool,
    /// Smallest constant lag, when known.
    pub min_constant_lag: Option<f64>,
    /// True for each reaction marked fast in the source.
    pub fast: Vec<bool>,
    /// Slots that hold `rateOf` values.
    pub rate_of_slots: Vec<usize>,
    /// Defining expressions of the initial state, re-evaluated for times before the start.
    pub initial_items: Vec<(usize, C)>,
    /// Rate per reaction with all derived symbols inlined, when inlining stayed within bounds.
    pub flat_rates: Vec<Expr>,
    /// Names of the state variables.
    pub state_names: Vec<String>,
}

const NODE_LIMIT: usize = 4_000_000;

/// Fast reactions run at this multiple of their law, a singular-perturbation stand-in for rapid
/// equilibrium: the residual of the equilibrium condition is of order `1 / FAST_SCALE`.
pub const FAST_SCALE: f64 = 1e8;

fn zero_rate_of(e: &Expr, is_local: &dyn Fn(&str) -> bool) -> Expr {
    match e {
        Expr::RateOf(s) if is_local(s) => Expr::Num(0.0),
        Expr::Num(_) | Expr::Sym(_) | Expr::Time | Expr::Avogadro | Expr::RateOf(_) => e.clone(),
        Expr::Apply(op, a) => Expr::Apply(*op, a.iter().map(|x| zero_rate_of(x, is_local)).collect()),
        Expr::Call(n, a) => Expr::Call(n.clone(), a.iter().map(|x| zero_rate_of(x, is_local)).collect()),
        Expr::Piecewise(p, o) => Expr::Piecewise(
            p.iter().map(|(v, c)| (zero_rate_of(v, is_local), zero_rate_of(c, is_local))).collect(),
            o.as_ref().map(|x| Box::new(zero_rate_of(x, is_local))),
        ),
        Expr::Delay(a, b) => Expr::Delay(Box::new(zero_rate_of(a, is_local)), Box::new(zero_rate_of(b, is_local))),
    }
}

fn inline_functions(e: &Expr, funcs: &HashMap<&str, &sbml::FunctionDefinition>) -> Result<Expr> {
    Ok(match e {
        Expr::Num(_) | Expr::Sym(_) | Expr::Time | Expr::Avogadro | Expr::RateOf(_) => e.clone(),
        Expr::Apply(op, args) => {
            Expr::Apply(*op, args.iter().map(|a| inline_functions(a, funcs)).collect::<Result<_>>()?)
        }
        Expr::Piecewise(pieces, otherwise) => Expr::Piecewise(
            pieces
                .iter()
                .map(|(v, c)| Ok((inline_functions(v, funcs)?, inline_functions(c, funcs)?)))
                .collect::<Result<_>>()?,
            match otherwise {
                Some(o) => Some(Box::new(inline_functions(o, funcs)?)),
                None => None,
            },
        ),
        Expr::Delay(a, b) => Expr::Delay(Box::new(inline_functions(a, funcs)?), Box::new(inline_functions(b, funcs)?)),
        Expr::Call(name, args) => {
            let f = funcs.get(name.as_str()).ok_or_else(|| Error::Kinetic(format!("unknown function `{name}`")))?;
            if f.args.len() != args.len() {
                return Err(Error::Kinetic(format!("wrong argument count for `{name}`")));
            }
            let actual: Vec<Expr> = args.iter().map(|a| inline_functions(a, funcs)).collect::<Result<_>>()?;
            let map: HashMap<&str, &Expr> = f.args.iter().map(String::as_str).zip(actual.iter()).collect();
            let body = inline_functions(&f.body, funcs)?;
            body.substitute(&|s| map.get(s).map(|e| (*e).clone()))
        }
    })
}

fn node_count(e: &Expr) -> usize {
    match e {
        Expr::Num(_) | Expr::Sym(_) | Expr::Time | Expr::Avogadro | Expr::RateOf(_) => 1,
        Expr::Apply(_, a) | Expr::Call(_, a) => 1 + a.iter().map(node_count).sum::<usize>(),
        Expr::Piecewise(p, o) => {
            1 + p.iter().map(|(v, c)| node_count(v) + node_count(c)).sum::<usize>()
                + o.as_ref().map_or(0, |e| node_count(e))
        }
        Expr::Delay(a, b) => 1 + node_count(a) + node_count(b),
    }
}

struct Builder<'a> {
    level: u32,
    slots: HashMap<String, usize>,
    names: Vec<String>,
    funcs: HashMap<&'a str, &'a sbml::FunctionDefinition>,
}

impl Builder<'_> {
    fn slot(&mut self, name: &str) -> usize {
        if let Some(&s) = self.slots.get(name) {
            return s;
        }
        let s = self.names.len();
        self.names.push(name.to_string());
        self.slots.insert(name.to_string(), s);
        s
    }

    fn inline(&self, e: &Expr) -> Result<Expr> {
        inline_functions(e, &self.funcs)
    }

    fn compile(&mut self, e: &Expr) -> Result<C> {
        compile_expr(e, &self.slots)
    }
}

/// Resolves the symbols of an inlined expression to slots.
pub fn compile_expr(e: &Expr, slots: &HashMap<String, usize>) -> Result<C> {
    Ok(match e {
        Expr::Num(v) => C::Num(*v),
        Expr::Time => C::Time,
        Expr::Avogadro => C::Num(AVOGADRO),
        Expr::Sym(s) => {
            let slot = *slots.get(s).ok_or_else(|| Error::Kinetic(format!("undefined symbol `{s}`")))?;
            C::Slot(slot)
        }
        Expr::RateOf(s) => {
            let name = format!("rateOf#{s}");
            let slot = *slots.get(&name).ok_or_else(|| Error::Kinetic(format!("no rate available for `{s}`")))?;
            C::Slot(slot)
        }
        Expr::Apply(op, args) => {
            let compiled: Vec<C> = args.iter().map(|a| compile_expr(a, slots)).collect::<Result<_>>()?;
            C::Ap(*op, compiled)
        }
        Expr::Piecewise(pieces, otherwise) => {
            let mut compiled = Vec::with_capacity(pieces.len());
            for (v, c) in pieces {
                compiled.push((compile_expr(v, slots)?, compile_expr(c, slots)?));
            }
            let otherwise = match otherwise {
                Some(o) => Some(Box::new(compile_expr(o, slots)?)),
                None => None,
            };
            C::Piece(compiled, otherwise)
        }
        Expr::Delay(a, b) => C::Delay(Box::new(compile_expr(a, slots)?), Box::new(compile_expr(b, slots)?)),
        Expr::Call(name, _) => return Err(Error::Kinetic(format!("function `{name}` was not inlined"))),
    })
}

fn collect_rate_of(e: &Expr, out: &mut Vec<String>) {
    match e {
        Expr::RateOf(s) => out.push(s.clone()),
        Expr::Num(_) | Expr::Sym(_) | Expr::Time | Expr::Avogadro => {}
        Expr::Apply(_, a) | Expr::Call(_, a) => a.iter().for_each(|x| collect_rate_of(x, out)),
        Expr::Piecewise(p, o) => {
            for (v, c) in p {
                collect_rate_of(v, out);
                collect_rate_of(c, out);
            }
            if let Some(x) = o {
                collect_rate_of(x, out);
            }
        }
        Expr::Delay(a, b) => {
            collect_rate_of(a, out);
            collect_rate_of(b, out);
        }
    }
}

fn constant_lags(e: &Expr, out: &mut Vec<f64>) {
    match e {
        Expr::Delay(a, b) => {
            if let Expr::Num(v) = **b {
                out.push(v);
            }
            constant_lags(a, out);
            constant_lags(b, out);
        }
        Expr::Num(_) | Expr::Sym(_) | Expr::Time | Expr::Avogadro | Expr::RateOf(_) => {}
        Expr::Apply(_, a) | Expr::Call(_, a) => a.iter().for_each(|x| constant_lags(x, out)),
        Expr::Piecewise(p, o) => {
            for (v, c) in p {
                constant_lags(v, out);
                constant_lags(c, out);
            }
            if let Some(x) = o {
                constant_lags(x, out);
            }
        }
    }
}

/// Finds which symbol each algebraic rule determines, by augmenting-path matching.
fn match_algebraic(rules: &[Expr], candidates: &HashSet<String>, prefer_not: &HashSet<&str>) -> Result<Vec<String>> {
    let options: Vec<Vec<String>> = rules
        .iter()
        .map(|r| {
            let mut syms = Vec::new();
            r.symbols(&mut syms);
            syms.sort();
            syms.dedup();
            syms.retain(|s| candidates.contains(s));
            // A symbol that events assign is a poor choice: the rule would undo the assignment.
            syms.sort_by_key(|s| prefer_not.contains(s.as_str()));
            syms
        })
        .collect();
    let mut owner: HashMap<String, usize> = HashMap::new();
    fn assign(
        rule: usize,
        options: &[Vec<String>],
        owner: &mut HashMap<String, usize>,
        seen: &mut HashSet<String>,
    ) -> bool {
        for sym in &options[rule] {
            if seen.insert(sym.clone()) {
                let holder = owner.get(sym).copied();
                let free = match holder {
                    None => true,
                    Some(other) => assign(other, options, owner, seen),
                };
                if free {
                    owner.insert(sym.clone(), rule);
                    return true;
                }
            }
        }
        false
    }
    for rule in 0..rules.len() {
        let mut seen = HashSet::new();
        if !assign(rule, &options, &mut owner, &mut seen) {
            return Err(Error::Kinetic("an algebraic rule has no free variable to determine".into()));
        }
    }
    let mut out = vec![String::new(); rules.len()];
    for (sym, rule) in owner {
        out[rule] = sym;
    }
    Ok(out)
}

impl KineticModel {
    /// Reads an SBML file and compiles it.
    pub fn read_sbml(path: impl AsRef<std::path::Path>) -> Result<Self> {
        Self::compile(&sbml::read::read_file(path)?)
    }

    /// Compiles SBML text.
    pub fn from_sbml_str(xml: &str) -> Result<Self> {
        Self::compile(&sbml::read::read_str(xml)?)
    }

    /// Compiles a parsed SBML document.
    pub fn compile(doc: &Document) -> Result<Self> {
        let src = &doc.model;
        let mut b = Builder {
            level: doc.level,
            slots: HashMap::new(),
            names: Vec::new(),
            funcs: src.functions.iter().map(|f| (f.id.as_str(), f)).collect(),
        };
        build(&mut b, src)
    }

    /// Overrides the initial value of a symbol and re-evaluates everything that depends on it.
    /// For a species name the value is an amount when written `S#a`, otherwise what math sees.
    pub fn set_value(&mut self, name: &str, value: f64) -> Result<()> {
        let slot = self.slot(name)?;
        self.initial_items.retain(|(s, _)| *s != slot);
        self.initial[slot] = value;
        // When the math value of an amount species is overridden, the amount follows.
        if let Some(sp) = self.species.iter().find(|s| s.view == slot)
            && matches!(sp.kind, SpeciesKind::Amount | SpeciesKind::Constant)
        {
            let size = self.initial[sp.compartment];
            let amount = if sp.has_only_substance_units { value } else { value * size };
            let amount_slot = sp.amount;
            self.initial_items.retain(|(s, _)| *s != amount_slot);
            self.initial[amount_slot] = amount;
        }
        for _ in 0..self.initial_items.len() + 3 {
            let mut changed = false;
            for (s, c) in &self.initial_items {
                let v = eval(c, &self.initial, 0.0, &NoDelay);
                let old = self.initial[*s];
                if v.to_bits() != old.to_bits() && !(v.is_nan() && old.is_nan()) {
                    self.initial[*s] = v;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        Ok(())
    }

    /// Number of state variables.
    pub fn n_states(&self) -> usize {
        self.state_slots.len()
    }

    /// The state at the initial time.
    pub fn initial_state(&self) -> Vec<f64> {
        self.state_slots.iter().map(|&s| self.initial[s]).collect()
    }

    /// Slot of a symbol.
    pub fn slot(&self, name: &str) -> Result<usize> {
        self.index.get(name).copied().ok_or_else(|| Error::Unknown { kind: "symbol", id: name.into() })
    }

    /// Fills `vals` from a state. `vals` must start as [`KineticModel::initial`] or a copy.
    pub fn eval_state(&self, t: f64, y: &[f64], vals: &mut [f64], ctx: &dyn DelayContext) {
        for (i, &s) in self.state_slots.iter().enumerate() {
            vals[s] = y[i];
        }
        for (slot, c) in &self.derived {
            vals[*slot] = eval(c, vals, t, ctx);
        }
    }

    /// Evaluates the right-hand side, assuming `vals` is current.
    pub fn rhs_from_vals(&self, t: f64, vals: &[f64], dy: &mut [f64], ctx: &dyn DelayContext) {
        dy.fill(0.0);
        for e in &self.stoich {
            let coef = match e.coef {
                Coef::Const(v) => v,
                Coef::Slot(s) => vals[s],
            };
            let conv = e.conversion.map_or(1.0, |s| vals[s]);
            dy[e.state] += e.sign * coef * conv * vals[self.reactions[e.reaction].rate];
        }
        for (i, c) in &self.rate_rules {
            dy[*i] = eval(c, vals, t, ctx);
        }
        for (i, c) in &self.algebraic {
            dy[*i] = eval(c, vals, t, ctx);
        }
    }
}

impl KineticModel {
    /// Right-hand side of the fast reactions alone, assuming `vals` is current.
    pub fn fast_rhs_from_vals(&self, vals: &[f64], dy: &mut [f64]) {
        dy.fill(0.0);
        for e in &self.stoich {
            if self.fast[e.reaction] {
                let coef = match e.coef {
                    Coef::Const(v) => v,
                    Coef::Slot(s) => vals[s],
                };
                let conv = e.conversion.map_or(1.0, |s| vals[s]);
                dy[e.state] += e.sign * coef * conv * vals[self.reactions[e.reaction].rate];
            }
        }
    }

    /// True if any reaction is fast.
    pub fn has_fast(&self) -> bool {
        self.fast.iter().any(|&f| f)
    }
}

/// Derived definition of a species' view and amount at run time and at the initial time.
struct SpeciesInit {
    view: Option<Expr>,
    amount: Option<Expr>,
}

fn build(b: &mut Builder<'_>, src: &SbmlModel) -> Result<KineticModel> {
    let local_substituted = |r: &sbml::Reaction| -> Result<Expr> {
        let Some(law) = &r.kinetic_law else { return Ok(Expr::Num(0.0)) };
        let locals: HashMap<&str, f64> =
            r.local_parameters.iter().map(|p| (p.id.as_str(), p.value.unwrap_or(f64::NAN))).collect();
        // A local parameter is constant, so `rateOf` of one is zero.
        let shadowed = zero_rate_of(law, &|s| locals.contains_key(s));
        let law = shadowed.substitute(&|s| locals.get(s).map(|v| Expr::Num(*v)));
        Ok(if r.fast { mul(law, Expr::Num(FAST_SCALE)) } else { law })
    };

    // Classify symbols.
    let mut assign_rule: HashMap<&str, Expr> = HashMap::new();
    let mut rate_rule: HashMap<&str, Expr> = HashMap::new();
    let mut algebraic_exprs: Vec<Expr> = Vec::new();
    for rule in &src.rules {
        match rule {
            Rule::Assignment { variable, math } => {
                assign_rule.insert(variable, b.inline(math)?);
            }
            Rule::Rate { variable, math } => {
                rate_rule.insert(variable, b.inline(math)?);
            }
            Rule::Algebraic { math } => algebraic_exprs.push(b.inline(math)?),
        }
    }
    let mut event_targets: HashSet<&str> = HashSet::new();
    for e in &src.events {
        for a in &e.assignments {
            event_targets.insert(&a.variable);
        }
    }
    let in_reaction: HashSet<&str> =
        src.reactions.iter().flat_map(|r| r.reactants.iter().chain(&r.products)).map(|s| s.species.as_str()).collect();

    // Candidates for algebraic variables: non-constant symbols no other rule or reaction determines.
    let mut candidates: HashSet<String> = HashSet::new();
    let free = |id: &str| !assign_rule.contains_key(id) && !rate_rule.contains_key(id);
    for c in &src.compartments {
        if !c.constant && free(&c.id) {
            candidates.insert(c.id.clone());
        }
    }
    for p in &src.parameters {
        if !p.constant && free(&p.id) {
            candidates.insert(p.id.clone());
        }
    }
    for s in &src.species {
        let driven = in_reaction.contains(s.id.as_str()) && !s.boundary_condition;
        if !s.constant && free(&s.id) && !driven {
            candidates.insert(s.id.clone());
        }
    }
    for r in &src.reactions {
        for sr in r.reactants.iter().chain(&r.products) {
            if let Some(id) = &sr.id
                && !sr.constant
                && free(id)
            {
                candidates.insert(id.clone());
            }
        }
    }
    let algebraic_vars = match_algebraic(&algebraic_exprs, &candidates, &event_targets)?;
    let algebraic_set: HashSet<&str> = algebraic_vars.iter().map(String::as_str).collect();

    // Slots of all named symbols.
    for c in &src.compartments {
        b.slot(&c.id);
    }
    for p in &src.parameters {
        b.slot(&p.id);
    }
    for s in &src.species {
        b.slot(&s.id);
        b.slot(&format!("{}#a", s.id));
    }
    for r in &src.reactions {
        b.slot(&r.id);
        for sr in r.reactants.iter().chain(&r.products) {
            if let Some(id) = &sr.id {
                b.slot(id);
            }
        }
    }

    // Every expression in the model, to discover rateOf and delay use.
    let mut laws: Vec<Expr> = Vec::with_capacity(src.reactions.len());
    for r in &src.reactions {
        laws.push(b.inline(&local_substituted(r)?)?);
    }
    let mut all_exprs: Vec<Expr> = laws.clone();
    all_exprs.extend(assign_rule.values().cloned());
    all_exprs.extend(rate_rule.values().cloned());
    all_exprs.extend(algebraic_exprs.iter().cloned());
    for ia in &src.initial_assignments {
        all_exprs.push(b.inline(&ia.math)?);
    }
    for e in &src.events {
        all_exprs.push(b.inline(&e.trigger)?);
        if let Some(d) = &e.delay {
            all_exprs.push(b.inline(d)?);
        }
        if let Some(p) = &e.priority {
            all_exprs.push(b.inline(p)?);
        }
        for a in &e.assignments {
            all_exprs.push(b.inline(&a.math)?);
        }
    }
    for r in &src.reactions {
        for sr in r.reactants.iter().chain(&r.products) {
            if let Some(math) = &sr.stoichiometry_math {
                all_exprs.push(b.inline(math)?);
            }
        }
    }
    let mut rate_of_targets: Vec<String> = Vec::new();
    for e in &all_exprs {
        collect_rate_of(e, &mut rate_of_targets);
    }
    let has_delay = all_exprs.iter().any(Expr::has_delay);
    let mut lags = Vec::new();
    for e in &all_exprs {
        constant_lags(e, &mut lags);
    }
    let min_constant_lag =
        lags.iter().copied().filter(|v| *v > 0.0).fold(None, |acc: Option<f64>, v| Some(acc.map_or(v, |a| a.min(v))));

    // Derived expressions and state variables.
    let mut derived_exprs: Vec<(String, Expr)> = Vec::new();
    let mut state_names: Vec<String> = Vec::new();
    let mut state_index: HashMap<String, usize> = HashMap::new();
    let mut add_state = |name: String, state_names: &mut Vec<String>| -> usize {
        let i = state_names.len();
        state_index.insert(name.clone(), i);
        state_names.push(name);
        i
    };

    // Non-species symbols.
    let mut value_state = |id: &str, state_names: &mut Vec<String>, derived: &mut Vec<(String, Expr)>| {
        if let Some(rule) = assign_rule.get(id) {
            derived.push((id.to_string(), rule.clone()));
        } else if rate_rule.contains_key(id) || algebraic_set.contains(id) || event_targets.contains(id) {
            add_state(id.to_string(), state_names);
        }
    };
    for c in &src.compartments {
        value_state(&c.id, &mut state_names, &mut derived_exprs);
    }
    for p in &src.parameters {
        value_state(&p.id, &mut state_names, &mut derived_exprs);
    }
    for r in &src.reactions {
        for sr in r.reactants.iter().chain(&r.products) {
            if let Some(id) = &sr.id {
                value_state(id, &mut state_names, &mut derived_exprs);
            }
        }
    }

    // Species.
    let mut species_infos = Vec::with_capacity(src.species.len());
    let mut species_init: Vec<SpeciesInit> = Vec::with_capacity(src.species.len());
    let ia_of: HashMap<&str, Expr> = {
        let mut map = HashMap::new();
        for ia in &src.initial_assignments {
            map.insert(ia.symbol.as_str(), b.inline(&ia.math)?);
        }
        map
    };
    for s in &src.species {
        let view_name = s.id.clone();
        let amount_name = format!("{}#a", s.id);
        let size = Expr::Sym(s.compartment.clone());
        let hosu = s.has_only_substance_units;
        let kind = if assign_rule.contains_key(s.id.as_str()) {
            SpeciesKind::Assigned
        } else if rate_rule.contains_key(s.id.as_str()) {
            SpeciesKind::RateRule
        } else if algebraic_set.contains(s.id.as_str()) {
            SpeciesKind::Algebraic
        } else if s.constant {
            SpeciesKind::Constant
        } else {
            SpeciesKind::Amount
        };
        let view_of_amount =
            if hosu { Expr::Sym(amount_name.clone()) } else { div(Expr::Sym(amount_name.clone()), size.clone()) };
        let amount_of_view =
            if hosu { Expr::Sym(view_name.clone()) } else { mul(Expr::Sym(view_name.clone()), size.clone()) };
        match kind {
            SpeciesKind::Amount => {
                add_state(amount_name.clone(), &mut state_names);
                derived_exprs.push((view_name.clone(), view_of_amount.clone()));
            }
            SpeciesKind::RateRule | SpeciesKind::Algebraic => {
                add_state(view_name.clone(), &mut state_names);
                derived_exprs.push((amount_name.clone(), amount_of_view.clone()));
            }
            SpeciesKind::Assigned => {
                derived_exprs.push((view_name.clone(), assign_rule[s.id.as_str()].clone()));
                derived_exprs.push((amount_name.clone(), amount_of_view.clone()));
            }
            SpeciesKind::Constant => {
                // The amount is what stays fixed when the compartment's size changes.
                derived_exprs.push((view_name.clone(), view_of_amount.clone()));
            }
        }
        // Initial definitions in order of precedence.
        let init = if let Some(rule) = assign_rule.get(s.id.as_str()) {
            SpeciesInit { view: Some(rule.clone()), amount: Some(amount_of_view.clone()) }
        } else if let Some(ia) = ia_of.get(s.id.as_str()) {
            SpeciesInit { view: Some(ia.clone()), amount: Some(amount_of_view.clone()) }
        } else if let Some(amount) = s.initial_amount {
            SpeciesInit { view: Some(view_of_amount.clone()), amount: Some(Expr::Num(amount)) }
        } else if let Some(conc) = s.initial_concentration {
            let view = if hosu { mul(Expr::Num(conc), size.clone()) } else { Expr::Num(conc) };
            SpeciesInit { view: Some(view), amount: Some(mul(Expr::Num(conc), size.clone())) }
        } else {
            SpeciesInit { view: None, amount: Some(amount_of_view.clone()) }
        };
        species_init.push(init);
        species_infos.push((s, kind));
    }

    // Reactions and stoichiometry.
    let mut reaction_infos = Vec::with_capacity(src.reactions.len());
    let mut stoich_exprs: Vec<(usize, usize, f64, Expr, Option<String>)> = Vec::new();
    let mut sm_counter = 0usize;
    let species_by_id: HashMap<&str, &sbml::Species> = src.species.iter().map(|s| (s.id.as_str(), s)).collect();
    let mut sr_initial: Vec<(String, Expr)> = Vec::new();
    for (k, r) in src.reactions.iter().enumerate() {
        derived_exprs.push((r.id.clone(), laws[k].clone()));
        reaction_infos.push(ReactionInfo { id: r.id.clone(), rate: b.slot(&r.id), law: laws[k].clone() });
        for (refs, sign) in [(&r.reactants, -1.0), (&r.products, 1.0)] {
            for sr in refs {
                let species = species_by_id
                    .get(sr.species.as_str())
                    .ok_or_else(|| Error::Unknown { kind: "species", id: sr.species.clone() })?;
                let coefficient = if let Some(math) = &sr.stoichiometry_math {
                    let name = format!("__stoich_math_{sm_counter}");
                    sm_counter += 1;
                    b.slot(&name);
                    derived_exprs.push((name.clone(), b.inline(math)?));
                    Expr::Sym(name)
                } else if let Some(id) = &sr.id {
                    sr_initial.push((
                        id.clone(),
                        Expr::Num(sr.stoichiometry.unwrap_or(if b.level < 3 { 1.0 } else { f64::NAN })),
                    ));
                    Expr::Sym(id.clone())
                } else {
                    Expr::Num(sr.stoichiometry.unwrap_or(1.0))
                };
                if species.boundary_condition || species.constant {
                    continue;
                }
                let name = format!("{}#a", species.id);
                let Some(&state) = state_index.get(&name) else { continue };
                let conversion = species.conversion_factor.clone().or_else(|| src.conversion_factor.clone());
                stoich_exprs.push((state, k, sign, coefficient, conversion));
            }
        }
    }

    // rateOf slots. Requested symbols get an expression built from the model's own definitions.
    let mut rate_defs: HashMap<String, Expr> = HashMap::new();
    let mut pending = rate_of_targets.clone();
    while let Some(target) = pending.pop() {
        if rate_defs.contains_key(&target) {
            continue;
        }
        let expr = rate_expression(&target, src, &assign_rule, &rate_rule, &stoich_exprs, &state_index, &mut pending);
        rate_defs.insert(target, expr);
    }
    for (name, expr) in &rate_defs {
        let slot_name = format!("rateOf#{name}");
        b.slot(&slot_name);
        derived_exprs.push((slot_name, expr.clone()));
    }

    // Conversion factor symbols are parameters already. Compile the state couplings.
    let state_slots: Vec<usize> = state_names.iter().map(|n| b.slot(n)).collect();
    let mut stoich = Vec::with_capacity(stoich_exprs.len());
    for (state, reaction, sign, coef, conversion) in &stoich_exprs {
        let coef = match coef {
            Expr::Num(v) => Coef::Const(*v),
            Expr::Sym(name) => Coef::Slot(b.slot(name)),
            _ => return Err(Error::Kinetic("unexpected stoichiometry expression".into())),
        };
        let conversion = match conversion {
            Some(c) => {
                Some(*b.slots.get(c).ok_or_else(|| Error::Unknown { kind: "conversion factor", id: c.clone() })?)
            }
            None => None,
        };
        stoich.push(StoichEntry { state: *state, reaction: *reaction, sign: *sign, coef, conversion });
    }

    // Order derived expressions by dependency.
    let derived_exprs = topo_sort(derived_exprs)?;
    let mut derived = Vec::with_capacity(derived_exprs.len());
    for (name, e) in &derived_exprs {
        derived.push((b.slot(name), b.compile(e)?));
    }

    let mut rate_rules = Vec::new();
    for (id, e) in &rate_rule {
        let state = state_position(&state_names, id)?;
        rate_rules.push((state, b.compile(e)?));
    }
    rate_rules.sort_by_key(|(i, _)| *i);
    let mut algebraic = Vec::new();
    for (k, e) in algebraic_exprs.iter().enumerate() {
        let state = state_position(&state_names, &algebraic_vars[k])?;
        algebraic.push((state, b.compile(e)?));
    }
    let mut mass = vec![1.0; state_names.len()];
    for (i, _) in &algebraic {
        mass[*i] = 0.0;
    }

    // Species records.
    let mut species = Vec::with_capacity(species_infos.len());
    for (s, kind) in &species_infos {
        species.push(SpeciesInfo {
            id: s.id.clone(),
            view: b.slots[s.id.as_str()],
            amount: b.slots[&format!("{}#a", s.id)],
            compartment: *b
                .slots
                .get(s.compartment.as_str())
                .ok_or_else(|| Error::Unknown { kind: "compartment", id: s.compartment.clone() })?,
            has_only_substance_units: s.has_only_substance_units,
            kind: *kind,
        });
    }

    // Events.
    let mut events = Vec::with_capacity(src.events.len());
    for e in &src.events {
        let mut assignments = Vec::with_capacity(e.assignments.len());
        for a in &e.assignments {
            let math = b.inline(&a.math)?;
            let compiled = b.compile(&math)?;
            let target = if let Some(sp) = species_by_id.get(a.variable.as_str()) {
                let amount_name = format!("{}#a", sp.id);
                if let Ok(state) = state_position(&state_names, &amount_name) {
                    Target::Amount {
                        state,
                        compartment: b.slots[sp.compartment.as_str()],
                        hosu: sp.has_only_substance_units,
                    }
                } else {
                    Target::State(state_position(&state_names, &a.variable)?)
                }
            } else {
                Target::State(state_position(&state_names, &a.variable)?)
            };
            assignments.push((target, compiled));
        }
        // Size changes come last: species values convert with the size before the event, and
        // rate-rule species then keep their amount across the jump in size.
        let is_size = |t: &Target| match t {
            Target::State(i) => src.compartments.iter().any(|c| c.id == state_names[*i]),
            Target::Amount { .. } => false,
        };
        assignments.sort_by_key(|(t, _)| is_size(t));
        let trigger = b.inline(&e.trigger)?;
        let delay = match &e.delay {
            Some(d) => Some(b.inline(d)?),
            None => None,
        };
        let priority = match &e.priority {
            Some(p) => Some(b.inline(p)?),
            None => None,
        };
        events.push(CompiledEvent {
            id: e.id.clone(),
            trigger: b.compile(&trigger)?,
            initial_value: e.initial_value,
            persistent: e.persistent,
            delay: match &delay {
                Some(d) => Some(b.compile(d)?),
                None => None,
            },
            priority: match &priority {
                Some(p) => Some(b.compile(p)?),
                None => None,
            },
            use_values_from_trigger_time: e.use_values_from_trigger_time,
            assignments,
        });
    }

    // Initial values: repeated evaluation of the defining expressions until stable.
    let mut init_items: Vec<(usize, C)> = Vec::new();
    let initial = {
        let mut plain: Vec<f64> = vec![f64::NAN; b.names.len()];
        for c in &src.compartments {
            let slot = b.slots[c.id.as_str()];
            if let Some(ia) = ia_of.get(c.id.as_str()) {
                init_items.push((slot, b.compile(ia)?));
            } else if !assign_rule.contains_key(c.id.as_str()) {
                plain[slot] = c.size.unwrap_or(if b.level < 3 { 1.0 } else { f64::NAN });
            }
        }
        for p in &src.parameters {
            let slot = b.slots[p.id.as_str()];
            if let Some(ia) = ia_of.get(p.id.as_str()) {
                init_items.push((slot, b.compile(ia)?));
            } else if !assign_rule.contains_key(p.id.as_str()) {
                plain[slot] = p.value.unwrap_or(f64::NAN);
            }
        }
        for (id, e) in &sr_initial {
            let slot = b.slots[id.as_str()];
            if let Some(ia) = ia_of.get(id.as_str()) {
                init_items.push((slot, b.compile(ia)?));
            } else if !assign_rule.contains_key(id.as_str())
                && let Expr::Num(v) = e
            {
                plain[slot] = *v;
            }
        }
        let mut species_overrides: HashSet<usize> = HashSet::new();
        for (k, s) in src.species.iter().enumerate() {
            let view_slot = b.slots[s.id.as_str()];
            let amount_slot = b.slots[&format!("{}#a", s.id)];
            species_overrides.insert(view_slot);
            species_overrides.insert(amount_slot);
            if let Some(e) = &species_init[k].view {
                init_items.push((view_slot, b.compile(e)?));
            }
            if let Some(e) = &species_init[k].amount {
                init_items.push((amount_slot, b.compile(e)?));
            }
        }
        // Run-time derived values participate unless a species definition above replaces them.
        for (slot, c) in &derived {
            if !species_overrides.contains(slot) {
                init_items.push((*slot, c.clone()));
            }
        }
        let mut vals = plain;
        let passes = init_items.len() + 3;
        let settle = |vals: &mut Vec<f64>, ctx: &dyn DelayContext| {
            for _ in 0..passes {
                let mut changed = false;
                for (slot, c) in &init_items {
                    let v = eval(c, vals, 0.0, ctx);
                    let old = vals[*slot];
                    if v.to_bits() != old.to_bits() && !(v.is_nan() && old.is_nan()) {
                        vals[*slot] = v;
                        changed = true;
                    }
                }
                if !changed {
                    break;
                }
            }
        };
        settle(&mut vals, &NoDelay);
        // A second round resolves `delay` in initial definitions, now that values exist to read.
        if has_delay {
            let base = vals.clone();
            let ctx = PreStart { items: &init_items, base: &base, rate_of_slots: &[] };
            settle(&mut vals, &ctx);
        }
        vals
    };

    let mut model = KineticModel {
        id: src.id.clone(),
        name: src.name.clone(),
        symbols: b.names.clone(),
        index: b.slots.clone(),
        state_slots,
        mass,
        species,
        reactions: reaction_infos,
        stoich,
        rate_rules,
        algebraic,
        events,
        derived,
        initial,
        flat_rhs: Vec::new(),
        jacobian: None,
        time_derivative: Vec::new(),
        has_delay,
        min_constant_lag,
        fast: src.reactions.iter().map(|r| r.fast).collect(),
        rate_of_slots: b.names.iter().enumerate().filter(|(_, n)| n.starts_with("rateOf#")).map(|(i, _)| i).collect(),
        initial_items: init_items,
        flat_rates: Vec::new(),
        state_names: state_names.clone(),
    };
    solve_algebraic_initial(&mut model);
    build_derivatives(
        &mut model,
        b,
        &derived_exprs,
        &stoich_exprs,
        &rate_rule,
        &algebraic_exprs,
        &algebraic_vars,
        &state_names,
    )?;
    Ok(model)
}

/// Values before the start of a simulation: the initial definitions (rules and initial assignments,
/// which may depend on time) re-evaluated at the earlier time. Rates are zero there.
pub struct PreStart<'a> {
    /// Defining expressions of the initial state.
    pub items: &'a [(usize, C)],
    /// Values to start from.
    pub base: &'a [f64],
    /// Slots of `rateOf` values.
    pub rate_of_slots: &'a [usize],
}

impl DelayContext for PreStart<'_> {
    fn value_at(&self, expr: &C, t: f64) -> f64 {
        let mut vals = self.base.to_vec();
        for &slot in self.rate_of_slots {
            vals[slot] = 0.0;
        }
        // Nested delays read the starting values, which bounds the depth of evaluation.
        for _ in 0..4 {
            for (slot, c) in self.items {
                if self.rate_of_slots.contains(slot) {
                    continue;
                }
                vals[*slot] = eval(c, &vals, t, &NoDelay);
            }
        }
        eval(expr, &vals, t, &NoDelay)
    }
}

fn state_position(state_names: &[String], name: &str) -> Result<usize> {
    state_names
        .iter()
        .position(|n| n == name)
        .ok_or_else(|| Error::Kinetic(format!("`{name}` is not a state variable")))
}

/// Builds the time derivative of a symbol from its definition.
fn rate_expression(
    target: &str,
    src: &SbmlModel,
    assign_rule: &HashMap<&str, Expr>,
    rate_rule: &HashMap<&str, Expr>,
    stoich: &[(usize, usize, f64, Expr, Option<String>)],
    state_index: &HashMap<String, usize>,
    pending: &mut Vec<String>,
) -> Expr {
    if let Some(e) = rate_rule.get(target) {
        return e.clone();
    }
    if let Some(e) = assign_rule.get(target) {
        // Chain rule through the symbols of the rule.
        let mut symbols = Vec::new();
        e.symbols(&mut symbols);
        symbols.sort();
        symbols.dedup();
        let mut total = diff::diff(e, Var::Time);
        for s in symbols {
            let d = diff::diff(e, Var::Sym(&s));
            if !matches!(d, Expr::Num(v) if v == 0.0) {
                pending.push(s.clone());
                total = add(total, mul(d, Expr::RateOf(s)));
            }
        }
        return total;
    }
    if let Some(sp) = src.species.iter().find(|s| s.id == target) {
        let name = format!("{}#a", sp.id);
        let Some(&state) = state_index.get(&name) else { return Expr::Num(0.0) };
        let mut amount_rate = Expr::Num(0.0);
        for (s, reaction, sign, coef, conversion) in stoich {
            if *s != state {
                continue;
            }
            let rate = Expr::Sym(src.reactions[*reaction].id.clone());
            let mut term = mul(mul(Expr::Num(*sign), coef.clone()), rate);
            if let Some(c) = conversion {
                term = mul(term, Expr::Sym(c.clone()));
            }
            amount_rate = add(amount_rate, term);
        }
        if sp.has_only_substance_units {
            return amount_rate;
        }
        let size = Expr::Sym(sp.compartment.clone());
        pending.push(sp.compartment.clone());
        return sub(
            div(amount_rate, size.clone()),
            mul(div(Expr::Sym(name), mul(size.clone(), size)), Expr::RateOf(sp.compartment.clone())),
        );
    }
    Expr::Num(0.0)
}

fn topo_sort(items: Vec<(String, Expr)>) -> Result<Vec<(String, Expr)>> {
    let position: HashMap<String, usize> = items.iter().enumerate().map(|(i, (n, _))| (n.clone(), i)).collect();
    let deps: Vec<Vec<usize>> = items
        .iter()
        .map(|(_, e)| {
            let mut syms = Vec::new();
            e.symbols(&mut syms);
            let mut rates = Vec::new();
            collect_rate_of(e, &mut rates);
            let mut out: Vec<usize> = syms.iter().filter_map(|s| position.get(s).copied()).collect();
            out.extend(rates.iter().filter_map(|s| position.get(&format!("rateOf#{s}")).copied()));
            out.sort_unstable();
            out.dedup();
            out
        })
        .collect();
    let mut state = vec![0u8; items.len()];
    let mut order = Vec::with_capacity(items.len());
    fn visit(
        i: usize,
        deps: &[Vec<usize>],
        state: &mut [u8],
        order: &mut Vec<usize>,
        names: &[(String, Expr)],
    ) -> Result<()> {
        match state[i] {
            2 => return Ok(()),
            1 => return Err(Error::Kinetic(format!("circular definition through `{}`", names[i].0))),
            _ => {}
        }
        state[i] = 1;
        for &d in &deps[i] {
            visit(d, deps, state, order, names)?;
        }
        state[i] = 2;
        order.push(i);
        Ok(())
    }
    for i in 0..items.len() {
        visit(i, &deps, &mut state, &mut order, &items)?;
    }
    let mut slots: Vec<Option<(String, Expr)>> = items.into_iter().map(Some).collect();
    Ok(order.into_iter().map(|i| slots[i].take().expect("visited once")).collect())
}

/// Solves algebraic rules for the variables they determine, by damped Newton on a numeric Jacobian.
fn solve_algebraic_initial(model: &mut KineticModel) {
    if model.algebraic.is_empty() {
        return;
    }
    let idx: Vec<usize> = model.algebraic.iter().map(|(i, _)| *i).collect();
    let mut y = model.initial_state();
    for &i in &idx {
        if !y[i].is_finite() {
            y[i] = 1.0;
        }
    }
    let n = idx.len();
    let mut vals = model.initial.clone();
    let residual = |y: &[f64], vals: &mut Vec<f64>, model: &KineticModel| -> Vec<f64> {
        model.eval_state(0.0, y, vals, &NoDelay);
        model.algebraic.iter().map(|(_, c)| eval(c, vals, 0.0, &NoDelay)).collect()
    };
    for _ in 0..50 {
        let r = residual(&y, &mut vals, model);
        let norm = r.iter().map(|v| v * v).sum::<f64>().sqrt();
        if !norm.is_finite() || norm < 1e-13 {
            break;
        }
        let mut jac = nalgebra::DMatrix::zeros(n, n);
        for (col, &i) in idx.iter().enumerate() {
            let save = y[i];
            let h = 1e-7 * save.abs().max(1e-3);
            y[i] = save + h;
            let plus = residual(&y, &mut vals, model);
            y[i] = save;
            for row in 0..n {
                jac[(row, col)] = (plus[row] - r[row]) / h;
            }
        }
        let rhs = nalgebra::DVector::from_iterator(n, r.iter().map(|v| -v));
        let Some(dx) = jac.lu().solve(&rhs) else { break };
        let mut scale = 1.0;
        for _ in 0..12 {
            let trial: Vec<f64> = {
                let mut t = y.clone();
                for (k, &i) in idx.iter().enumerate() {
                    t[i] += scale * dx[k];
                }
                t
            };
            let rt = residual(&trial, &mut vals, model);
            let nt = rt.iter().map(|v| v * v).sum::<f64>().sqrt();
            if nt.is_finite() && nt < norm {
                y = trial;
                break;
            }
            scale *= 0.5;
        }
    }
    model.eval_state(0.0, &y, &mut vals, &NoDelay);
    model.initial = vals;
}

#[allow(clippy::too_many_arguments)]
fn build_derivatives(
    model: &mut KineticModel,
    b: &mut Builder<'_>,
    derived_exprs: &[(String, Expr)],
    stoich_exprs: &[(usize, usize, f64, Expr, Option<String>)],
    rate_rule: &HashMap<&str, Expr>,
    algebraic_exprs: &[Expr],
    algebraic_vars: &[String],
    state_names: &[String],
) -> Result<()> {
    // Inline the derived definitions into each other, in evaluation order.
    let mut flat: HashMap<String, Expr> = HashMap::new();
    let mut budget = NODE_LIMIT;
    let expand_rates = |e: &Expr, flat: &HashMap<String, Expr>| -> Expr { expand_rate_of(e, flat) };
    for (name, e) in derived_exprs {
        let inlined = expand_rates(e, &flat).substitute(&|s| flat.get(s).cloned());
        let count = node_count(&inlined);
        if count > budget {
            return Ok(());
        }
        budget -= count;
        flat.insert(name.clone(), inlined);
    }
    let mut rhs: Vec<Expr> = vec![Expr::Num(0.0); state_names.len()];
    for (state, reaction, sign, coef, conversion) in stoich_exprs {
        let rate = flat.get(&model.reactions[*reaction].id).cloned().unwrap_or(Expr::Num(0.0));
        let coef_flat = coef.substitute(&|s| flat.get(s).cloned());
        let mut term = mul(mul(Expr::Num(*sign), coef_flat), rate);
        if let Some(c) = conversion {
            let cf = Expr::Sym(c.clone()).substitute(&|s| flat.get(s).cloned());
            term = mul(term, cf);
        }
        let old = std::mem::replace(&mut rhs[*state], Expr::Num(0.0));
        rhs[*state] = add(old, term);
    }
    for (id, e) in rate_rule {
        if let Some(i) = state_names.iter().position(|n| n == id) {
            rhs[i] = expand_rates(e, &flat).substitute(&|s| flat.get(s).cloned());
        }
    }
    for (k, e) in algebraic_exprs.iter().enumerate() {
        if let Some(i) = state_names.iter().position(|n| *n == algebraic_vars[k]) {
            rhs[i] = expand_rates(e, &flat).substitute(&|s| flat.get(s).cloned());
        }
    }
    let state_set: HashMap<&str, usize> = state_names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut jacobian = Vec::new();
    let mut time_derivative = Vec::new();
    for (i, f) in rhs.iter().enumerate() {
        let mut syms = Vec::new();
        f.symbols(&mut syms);
        syms.sort();
        syms.dedup();
        for s in syms {
            if let Some(&j) = state_set.get(s.as_str()) {
                let d = diff::diff(f, Var::Sym(&s));
                if !matches!(d, Expr::Num(v) if v == 0.0) {
                    if node_count(&d) > budget {
                        return Ok(());
                    }
                    budget -= node_count(&d);
                    jacobian.push((i, j, b.compile(&d)?));
                }
            }
        }
        let dt = diff::diff(f, Var::Time);
        if !matches!(dt, Expr::Num(v) if v == 0.0) {
            time_derivative.push((i, b.compile(&dt)?));
        }
    }
    let _ = neg;
    model.flat_rates = model.reactions.iter().map(|r| flat.get(&r.id).cloned().unwrap_or(Expr::Num(0.0))).collect();
    model.flat_rhs = rhs;
    model.jacobian = Some(jacobian);
    model.time_derivative = time_derivative;
    Ok(())
}

/// Replaces `rateOf(x)` by the flattened definition of its slot, when one exists.
fn expand_rate_of(e: &Expr, flat: &HashMap<String, Expr>) -> Expr {
    match e {
        Expr::RateOf(s) => flat.get(&format!("rateOf#{s}")).cloned().unwrap_or(Expr::Num(0.0)),
        Expr::Num(_) | Expr::Sym(_) | Expr::Time | Expr::Avogadro => e.clone(),
        Expr::Apply(op, a) => Expr::Apply(*op, a.iter().map(|x| expand_rate_of(x, flat)).collect()),
        Expr::Call(n, a) => Expr::Call(n.clone(), a.iter().map(|x| expand_rate_of(x, flat)).collect()),
        Expr::Piecewise(p, o) => Expr::Piecewise(
            p.iter().map(|(v, c)| (expand_rate_of(v, flat), expand_rate_of(c, flat))).collect(),
            o.as_ref().map(|x| Box::new(expand_rate_of(x, flat))),
        ),
        Expr::Delay(a, b) => Expr::Delay(Box::new(expand_rate_of(a, flat)), Box::new(expand_rate_of(b, flat))),
    }
}
