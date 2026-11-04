//! Mathematical expressions of SBML kinetic models, with evaluation and symbolic differentiation.

use serde::{Deserialize, Serialize};

/// Built-in operators of SBML MathML.
#[allow(missing_docs)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Hash)]
pub enum Op {
    Add,
    Sub,
    Mul,
    Div,
    Neg,
    Pow,
    Root,
    Abs,
    Exp,
    Ln,
    Log,
    Floor,
    Ceil,
    Factorial,
    Eq,
    Neq,
    Gt,
    Lt,
    Geq,
    Leq,
    And,
    Or,
    Xor,
    Not,
    Implies,
    Sin,
    Cos,
    Tan,
    Sec,
    Csc,
    Cot,
    Sinh,
    Cosh,
    Tanh,
    Sech,
    Csch,
    Coth,
    Asin,
    Acos,
    Atan,
    Asec,
    Acsc,
    Acot,
    Asinh,
    Acosh,
    Atanh,
    Asech,
    Acsch,
    Acoth,
    Max,
    Min,
    Quotient,
    Rem,
}

/// An expression tree. Booleans evaluate to 1 and 0.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Expr {
    /// Numeric constant.
    Num(f64),
    /// Named symbol (species, compartment, parameter, species reference or local parameter).
    Sym(String),
    /// Simulation time.
    Time,
    /// Avogadro's constant.
    Avogadro,
    /// Built-in operator application.
    Apply(Op, Vec<Expr>),
    /// `(value, condition)` pieces followed by the optional otherwise branch.
    Piecewise(Vec<(Expr, Expr)>, Option<Box<Expr>>),
    /// Call of a function definition.
    Call(String, Vec<Expr>),
    /// `delay(expression, lag)`.
    Delay(Box<Expr>, Box<Expr>),
    /// `rateOf(symbol)`.
    RateOf(String),
}

impl Expr {
    /// Binary application.
    pub fn bin(op: Op, a: Expr, b: Expr) -> Expr {
        Expr::Apply(op, vec![a, b])
    }

    /// True if the expression refers to the symbol `name`.
    pub fn mentions(&self, name: &str) -> bool {
        match self {
            Expr::Sym(s) | Expr::RateOf(s) => s == name,
            Expr::Num(_) | Expr::Time | Expr::Avogadro => false,
            Expr::Apply(_, a) | Expr::Call(_, a) => a.iter().any(|e| e.mentions(name)),
            Expr::Piecewise(p, o) => {
                p.iter().any(|(v, c)| v.mentions(name) || c.mentions(name))
                    || o.as_ref().is_some_and(|e| e.mentions(name))
            }
            Expr::Delay(a, b) => a.mentions(name) || b.mentions(name),
        }
    }

    /// Collects every symbol the expression refers to.
    pub fn symbols(&self, out: &mut Vec<String>) {
        match self {
            Expr::Sym(s) | Expr::RateOf(s) => out.push(s.clone()),
            Expr::Num(_) | Expr::Time | Expr::Avogadro => {}
            Expr::Apply(_, a) | Expr::Call(_, a) => a.iter().for_each(|e| e.symbols(out)),
            Expr::Piecewise(p, o) => {
                for (v, c) in p {
                    v.symbols(out);
                    c.symbols(out);
                }
                if let Some(e) = o {
                    e.symbols(out);
                }
            }
            Expr::Delay(a, b) => {
                a.symbols(out);
                b.symbols(out);
            }
        }
    }

    /// True if the expression contains a `delay` call.
    pub fn has_delay(&self) -> bool {
        match self {
            Expr::Delay(..) => true,
            Expr::Num(_) | Expr::Sym(_) | Expr::Time | Expr::Avogadro | Expr::RateOf(_) => false,
            Expr::Apply(_, a) | Expr::Call(_, a) => a.iter().any(Expr::has_delay),
            Expr::Piecewise(p, o) => {
                p.iter().any(|(v, c)| v.has_delay() || c.has_delay()) || o.as_ref().is_some_and(|e| e.has_delay())
            }
        }
    }

    /// Replaces symbols through `f`; symbols for which `f` returns `None` stay.
    pub fn substitute(&self, f: &dyn Fn(&str) -> Option<Expr>) -> Expr {
        match self {
            Expr::Sym(s) => f(s).unwrap_or_else(|| self.clone()),
            Expr::Num(_) | Expr::Time | Expr::Avogadro | Expr::RateOf(_) => self.clone(),
            Expr::Apply(op, a) => Expr::Apply(*op, a.iter().map(|e| e.substitute(f)).collect()),
            Expr::Call(n, a) => Expr::Call(n.clone(), a.iter().map(|e| e.substitute(f)).collect()),
            Expr::Piecewise(p, o) => Expr::Piecewise(
                p.iter().map(|(v, c)| (v.substitute(f), c.substitute(f))).collect(),
                o.as_ref().map(|e| Box::new(e.substitute(f))),
            ),
            Expr::Delay(a, b) => Expr::Delay(Box::new(a.substitute(f)), Box::new(b.substitute(f))),
        }
    }
}

/// Applies a built-in operator to already evaluated arguments.
pub fn apply_op(op: Op, a: &[f64]) -> f64 {
    let b = |v: bool| if v { 1.0 } else { 0.0 };
    let t = |v: f64| v != 0.0;
    match op {
        Op::Add => a.iter().sum(),
        Op::Mul => a.iter().product(),
        Op::Sub => {
            if a.len() == 1 {
                -a[0]
            } else {
                a[0] - a[1]
            }
        }
        Op::Neg => -a[0],
        Op::Div => a[0] / a[1],
        Op::Pow => pow(a[0], a[1]),
        Op::Root => {
            // root(degree, radicand)
            if a.len() == 1 {
                a[0].sqrt()
            } else if a[0] == 2.0 {
                a[1].sqrt()
            } else if a[0] == 3.0 {
                a[1].cbrt()
            } else if a[1] < 0.0 && a[0].fract() == 0.0 && (a[0] as i64) % 2 != 0 {
                -pow(-a[1], 1.0 / a[0])
            } else {
                pow(a[1], 1.0 / a[0])
            }
        }
        Op::Abs => a[0].abs(),
        Op::Exp => a[0].exp(),
        Op::Ln => a[0].ln(),
        Op::Log => {
            // log(base, value)
            if a.len() == 1 {
                a[0].log10()
            } else if a[0] == 10.0 {
                a[1].log10()
            } else if a[0] == 2.0 {
                a[1].log2()
            } else {
                a[1].ln() / a[0].ln()
            }
        }
        Op::Floor => a[0].floor(),
        Op::Ceil => a[0].ceil(),
        Op::Factorial => factorial(a[0]),
        Op::Eq => b(a.windows(2).all(|w| w[0] == w[1])),
        Op::Neq => b(a[0] != a[1]),
        Op::Gt => b(a.windows(2).all(|w| w[0] > w[1])),
        Op::Lt => b(a.windows(2).all(|w| w[0] < w[1])),
        Op::Geq => b(a.windows(2).all(|w| w[0] >= w[1])),
        Op::Leq => b(a.windows(2).all(|w| w[0] <= w[1])),
        Op::And => b(a.iter().all(|&v| t(v))),
        Op::Or => b(a.iter().any(|&v| t(v))),
        Op::Xor => b(a.iter().filter(|&&v| t(v)).count() % 2 == 1),
        Op::Not => b(!t(a[0])),
        Op::Implies => b(!t(a[0]) || t(a[1])),
        Op::Sin => a[0].sin(),
        Op::Cos => a[0].cos(),
        Op::Tan => a[0].tan(),
        Op::Sec => 1.0 / a[0].cos(),
        Op::Csc => 1.0 / a[0].sin(),
        Op::Cot => 1.0 / a[0].tan(),
        Op::Sinh => a[0].sinh(),
        Op::Cosh => a[0].cosh(),
        Op::Tanh => a[0].tanh(),
        Op::Sech => 1.0 / a[0].cosh(),
        Op::Csch => 1.0 / a[0].sinh(),
        Op::Coth => 1.0 / a[0].tanh(),
        Op::Asin => a[0].asin(),
        Op::Acos => a[0].acos(),
        Op::Atan => a[0].atan(),
        Op::Asec => (1.0 / a[0]).acos(),
        Op::Acsc => (1.0 / a[0]).asin(),
        Op::Acot => (1.0 / a[0]).atan(),
        Op::Asinh => a[0].asinh(),
        Op::Acosh => a[0].acosh(),
        Op::Atanh => a[0].atanh(),
        Op::Asech => (1.0 / a[0]).acosh(),
        Op::Acsch => (1.0 / a[0]).asinh(),
        Op::Acoth => (1.0 / a[0]).atanh(),
        Op::Max => a.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        Op::Min => a.iter().copied().fold(f64::INFINITY, f64::min),
        Op::Quotient => (a[0] / a[1]).trunc(),
        Op::Rem => a[0] % a[1],
    }
}

fn pow(x: f64, y: f64) -> f64 {
    if y == 2.0 {
        x * x
    } else if y == 1.0 {
        x
    } else if y.fract() == 0.0 && y.abs() < 64.0 {
        x.powi(y as i32)
    } else {
        x.powf(y)
    }
}

fn factorial(x: f64) -> f64 {
    if x < 0.0 || x.fract() != 0.0 {
        return f64::NAN;
    }
    if x > 170.0 {
        return f64::INFINITY;
    }
    (1..=(x as u64)).map(|v| v as f64).product()
}
