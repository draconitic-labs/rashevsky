//! Symbolic differentiation of expressions, with constant folding at construction.

use super::expr::{Expr, Op};

/// Variable of differentiation.
#[derive(Clone, Copy, Debug)]
pub enum Var<'a> {
    /// A named symbol.
    Sym(&'a str),
    /// Simulation time.
    Time,
}

fn num(v: f64) -> Expr {
    Expr::Num(v)
}

fn is_zero(e: &Expr) -> bool {
    matches!(e, Expr::Num(v) if *v == 0.0)
}

fn is_one(e: &Expr) -> bool {
    matches!(e, Expr::Num(v) if *v == 1.0)
}

/// `a + b` with folding.
pub fn add(a: Expr, b: Expr) -> Expr {
    match (&a, &b) {
        (Expr::Num(x), Expr::Num(y)) => num(x + y),
        _ if is_zero(&a) => b,
        _ if is_zero(&b) => a,
        _ => Expr::bin(Op::Add, a, b),
    }
}

/// `a - b` with folding.
pub fn sub(a: Expr, b: Expr) -> Expr {
    match (&a, &b) {
        (Expr::Num(x), Expr::Num(y)) => num(x - y),
        _ if is_zero(&b) => a,
        _ if is_zero(&a) => neg(b),
        _ => Expr::bin(Op::Sub, a, b),
    }
}

/// `-a` with folding.
pub fn neg(a: Expr) -> Expr {
    match a {
        Expr::Num(x) => num(-x),
        Expr::Apply(Op::Neg, mut v) => v.remove(0),
        other => Expr::Apply(Op::Neg, vec![other]),
    }
}

/// `a * b` with folding.
pub fn mul(a: Expr, b: Expr) -> Expr {
    match (&a, &b) {
        (Expr::Num(x), Expr::Num(y)) => num(x * y),
        _ if is_zero(&a) || is_zero(&b) => num(0.0),
        _ if is_one(&a) => b,
        _ if is_one(&b) => a,
        _ => Expr::bin(Op::Mul, a, b),
    }
}

/// `a / b` with folding.
pub fn div(a: Expr, b: Expr) -> Expr {
    match (&a, &b) {
        (Expr::Num(x), Expr::Num(y)) => num(x / y),
        _ if is_zero(&a) => num(0.0),
        _ if is_one(&b) => a,
        _ => Expr::bin(Op::Div, a, b),
    }
}

/// `a ^ b` with folding.
pub fn pow(a: Expr, b: Expr) -> Expr {
    match (&a, &b) {
        (Expr::Num(x), Expr::Num(y)) => num(x.powf(*y)),
        _ if is_zero(&b) => num(1.0),
        _ if is_one(&b) => a,
        _ => Expr::bin(Op::Pow, a, b),
    }
}

fn ap1(op: Op, a: Expr) -> Expr {
    Expr::Apply(op, vec![a])
}

fn sign_of(a: &Expr) -> Expr {
    Expr::Piecewise(
        vec![(num(1.0), Expr::bin(Op::Gt, a.clone(), num(0.0))), (num(-1.0), Expr::bin(Op::Lt, a.clone(), num(0.0)))],
        Some(Box::new(num(0.0))),
    )
}

/// Derivative of `e` with respect to `x`.
///
/// Piecewise functions differentiate branch by branch. Step-like operators (`floor`, comparisons,
/// Boolean operators, `factorial`, `quotient`, `rem`) have zero derivative almost everywhere.
/// `delay` is treated as independent of the current state, and `rateOf` must be expanded before
/// differentiation; both differentiate to zero here.
pub fn diff(e: &Expr, x: Var<'_>) -> Expr {
    match e {
        Expr::Num(_) | Expr::Avogadro | Expr::Delay(..) | Expr::RateOf(_) | Expr::Call(..) => num(0.0),
        Expr::Time => num(if matches!(x, Var::Time) { 1.0 } else { 0.0 }),
        Expr::Sym(s) => num(if matches!(x, Var::Sym(v) if v == s) { 1.0 } else { 0.0 }),
        Expr::Piecewise(pieces, otherwise) => Expr::Piecewise(
            pieces.iter().map(|(v, c)| (diff(v, x), c.clone())).collect(),
            otherwise.as_ref().map(|o| Box::new(diff(o, x))),
        ),
        Expr::Apply(op, args) => diff_apply(*op, args, x),
    }
}

fn diff_apply(op: Op, a: &[Expr], x: Var<'_>) -> Expr {
    let d = |k: usize| diff(&a[k], x);
    let chain = |outer: Expr| mul(outer, diff(&a[0], x));
    let u = || a[0].clone();
    match op {
        Op::Add => a.iter().map(|t| diff(t, x)).fold(num(0.0), add),
        Op::Sub => {
            if a.len() == 1 {
                neg(d(0))
            } else {
                sub(d(0), d(1))
            }
        }
        Op::Neg => neg(d(0)),
        Op::Mul => {
            // Product rule over n factors.
            let mut total = num(0.0);
            for k in 0..a.len() {
                let dk = d(k);
                if is_zero(&dk) {
                    continue;
                }
                let mut term = dk;
                for (m, factor) in a.iter().enumerate() {
                    if m != k {
                        term = mul(term, factor.clone());
                    }
                }
                total = add(total, term);
            }
            total
        }
        Op::Div => {
            let (da, db) = (d(0), d(1));
            if is_zero(&db) {
                div(da, a[1].clone())
            } else {
                div(
                    sub(mul(da, a[1].clone()), mul(a[0].clone(), db)),
                    mul(a[1].clone(), a[1].clone()),
                )
            }
        }
        Op::Pow => power_rule(&a[0], &a[1], x),
        Op::Root => {
            if a.len() == 1 {
                power_rule(&a[0], &num(0.5), x)
            } else {
                power_rule(&a[1], &div(num(1.0), a[0].clone()), x)
            }
        }
        Op::Abs => chain(sign_of(&a[0])),
        Op::Exp => chain(ap1(Op::Exp, u())),
        Op::Ln => div(d(0), u()),
        Op::Log => {
            // log(base, value) with constant or variable base.
            let (base, value) = if a.len() == 1 { (num(10.0), a[0].clone()) } else { (a[0].clone(), a[1].clone()) };
            let ratio = div(ap1(Op::Ln, value), ap1(Op::Ln, base));
            diff(&ratio, x)
        }
        Op::Sin => chain(ap1(Op::Cos, u())),
        Op::Cos => chain(neg(ap1(Op::Sin, u()))),
        Op::Tan => chain(pow(ap1(Op::Sec, u()), num(2.0))),
        Op::Sec => chain(mul(ap1(Op::Sec, u()), ap1(Op::Tan, u()))),
        Op::Csc => chain(neg(mul(ap1(Op::Csc, u()), ap1(Op::Cot, u())))),
        Op::Cot => chain(neg(pow(ap1(Op::Csc, u()), num(2.0)))),
        Op::Sinh => chain(ap1(Op::Cosh, u())),
        Op::Cosh => chain(ap1(Op::Sinh, u())),
        Op::Tanh => chain(pow(ap1(Op::Sech, u()), num(2.0))),
        Op::Sech => chain(neg(mul(ap1(Op::Sech, u()), ap1(Op::Tanh, u())))),
        Op::Csch => chain(neg(mul(ap1(Op::Csch, u()), ap1(Op::Coth, u())))),
        Op::Coth => chain(neg(pow(ap1(Op::Csch, u()), num(2.0)))),
        Op::Asin => chain(div(num(1.0), root2(sub(num(1.0), square(u()))))),
        Op::Acos => chain(neg(div(num(1.0), root2(sub(num(1.0), square(u())))))),
        Op::Atan => chain(div(num(1.0), add(num(1.0), square(u())))),
        Op::Acot => chain(neg(div(num(1.0), add(num(1.0), square(u()))))),
        Op::Asec => chain(div(num(1.0), mul(ap1(Op::Abs, u()), root2(sub(square(u()), num(1.0)))))),
        Op::Acsc => chain(neg(div(num(1.0), mul(ap1(Op::Abs, u()), root2(sub(square(u()), num(1.0))))))),
        Op::Asinh => chain(div(num(1.0), root2(add(square(u()), num(1.0))))),
        Op::Acosh => chain(div(num(1.0), root2(sub(square(u()), num(1.0))))),
        Op::Atanh | Op::Acoth => chain(div(num(1.0), sub(num(1.0), square(u())))),
        Op::Asech => chain(neg(div(num(1.0), mul(u(), root2(sub(num(1.0), square(u()))))))),
        Op::Acsch => chain(neg(div(num(1.0), mul(ap1(Op::Abs, u()), root2(add(num(1.0), square(u()))))))),
        Op::Max | Op::Min => {
            // Reduce to a chain of two-way branches.
            let cmp = if op == Op::Max { Op::Geq } else { Op::Leq };
            let mut best = a[a.len() - 1].clone();
            let mut derivative = diff(&best, x);
            for item in a[..a.len() - 1].iter().rev() {
                let cond = Expr::bin(cmp, item.clone(), best.clone());
                derivative = Expr::Piecewise(vec![(diff(item, x), cond.clone())], Some(Box::new(derivative)));
                best = Expr::Piecewise(vec![(item.clone(), cond)], Some(Box::new(best)));
            }
            derivative
        }
        Op::Floor | Op::Ceil | Op::Factorial | Op::Eq | Op::Neq | Op::Gt | Op::Lt | Op::Geq | Op::Leq
        | Op::And | Op::Or | Op::Xor | Op::Not | Op::Implies | Op::Quotient | Op::Rem => num(0.0),
    }
}

fn square(e: Expr) -> Expr {
    pow(e, num(2.0))
}

fn root2(e: Expr) -> Expr {
    Expr::Apply(Op::Root, vec![num(2.0), e])
}

fn power_rule(base: &Expr, exponent: &Expr, x: Var<'_>) -> Expr {
    let db = diff(base, x);
    let de = diff(exponent, x);
    if is_zero(&de) {
        if is_zero(&db) {
            return num(0.0);
        }
        // d(b^c) = c b^(c-1) db
        let reduced = match exponent {
            Expr::Num(c) => num(c - 1.0),
            other => sub(other.clone(), num(1.0)),
        };
        return mul(mul(exponent.clone(), pow(base.clone(), reduced)), db);
    }
    // d(b^e) = b^e (de ln b + e db / b)
    let whole = pow(base.clone(), exponent.clone());
    let left = mul(de, ap1(Op::Ln, base.clone()));
    let right = if is_zero(&db) { num(0.0) } else { div(mul(exponent.clone(), db), base.clone()) };
    mul(whole, add(left, right))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kinetic::expr::apply_op;
    use proptest::prelude::*;

    fn eval(e: &Expr, x: f64) -> f64 {
        match e {
            Expr::Num(v) => *v,
            Expr::Sym(_) => x,
            Expr::Time => 0.0,
            Expr::Apply(op, args) => {
                let v: Vec<f64> = args.iter().map(|a| eval(a, x)).collect();
                apply_op(*op, &v)
            }
            Expr::Piecewise(pieces, o) => {
                for (v, c) in pieces {
                    if eval(c, x) != 0.0 {
                        return eval(v, x);
                    }
                }
                o.as_ref().map_or(f64::NAN, |e| eval(e, x))
            }
            _ => f64::NAN,
        }
    }

    fn sym() -> Expr {
        Expr::Sym("x".into())
    }

    fn check_all(x: f64) {
        let s = sym;
        let ops = [
            Op::Sin, Op::Cos, Op::Tan, Op::Sec, Op::Csc, Op::Cot, Op::Sinh, Op::Cosh, Op::Tanh, Op::Sech,
            Op::Csch, Op::Coth, Op::Atan, Op::Acot, Op::Asinh, Op::Exp, Op::Ln, Op::Abs,
        ];
        let mut cases: Vec<Expr> = ops.iter().map(|o| Expr::Apply(*o, vec![s()])).collect();
        cases.push(Expr::bin(Op::Pow, s(), Expr::Num(3.0)));
        cases.push(Expr::bin(Op::Pow, Expr::Num(2.0), s()));
        cases.push(Expr::bin(Op::Pow, s(), s()));
        cases.push(Expr::bin(Op::Div, Expr::Apply(Op::Sin, vec![s()]), Expr::bin(Op::Add, s(), Expr::Num(2.0))));
        cases.push(Expr::Apply(Op::Mul, vec![s(), s(), Expr::Apply(Op::Exp, vec![s()])]));
        cases.push(Expr::Apply(Op::Log, vec![Expr::Num(3.0), s()]));
        cases.push(Expr::Apply(Op::Root, vec![Expr::Num(3.0), s()]));
        cases.push(Expr::Apply(Op::Max, vec![s(), Expr::Num(1.5), Expr::bin(Op::Mul, Expr::Num(2.0), s())]));
        for e in cases {
            let d = diff(&e, Var::Sym("x"));
            let h = 1e-6 * x.abs().max(1.0);
            let numeric = (eval(&e, x + h) - eval(&e, x - h)) / (2.0 * h);
            let analytic = eval(&d, x);
            assert!(
                (numeric - analytic).abs() < 1e-5 * (1.0 + numeric.abs()),
                "{e:?} at {x}: numeric {numeric}, analytic {analytic}"
            );
        }
    }

    #[test]
    fn matches_central_differences_at_fixed_points() {
        for x in [0.4, 0.73, 1.3, 2.7] {
            check_all(x);
        }
    }

    proptest! {
        #[test]
        fn matches_central_differences(x in 0.2f64..3.0) {
            // Avoid the poles of the reciprocal trigonometric functions.
            prop_assume!((x - std::f64::consts::FRAC_PI_2).abs() > 0.1);
            prop_assume!((x - std::f64::consts::PI).abs() > 0.1);
            check_all(x);
        }
    }

    #[test]
    fn inverse_functions_inside_their_domains() {
        let s = sym;
        for (op, x) in [(Op::Asin, 0.3), (Op::Acos, 0.3), (Op::Atanh, 0.3), (Op::Asech, 0.4), (Op::Acsch, 0.8)] {
            let e = Expr::Apply(op, vec![s()]);
            let d = diff(&e, Var::Sym("x"));
            let h = 1e-6;
            let numeric = (eval(&e, x + h) - eval(&e, x - h)) / (2.0 * h);
            assert!((numeric - eval(&d, x)).abs() < 1e-6, "{op:?}");
        }
        for (op, x) in [(Op::Asec, 1.7), (Op::Acsc, 1.7), (Op::Acosh, 1.7), (Op::Acoth, 1.7)] {
            let e = Expr::Apply(op, vec![s()]);
            let d = diff(&e, Var::Sym("x"));
            let h = 1e-6;
            let numeric = (eval(&e, x + h) - eval(&e, x - h)) / (2.0 * h);
            assert!((numeric - eval(&d, x)).abs() < 1e-6, "{op:?}");
        }
    }
}
