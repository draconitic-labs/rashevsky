//! MathML (the SBML subset) to and from [`Expr`].

use super::dom::Node;
use crate::error::{Error, Result};
use crate::kinetic::expr::{Expr, Op};

fn bad(message: impl Into<String>) -> Error {
    Error::Sbml(format!("mathml: {}", message.into()))
}

/// Reads the expression inside a `<math>` element.
pub fn read_math(math: &Node) -> Result<Expr> {
    let inner = math.children.first().ok_or_else(|| bad("empty <math>"))?;
    read_expr(strip_prefix(inner))
}

/// Reads a `<math>` that holds a `lambda`, returning argument names and body.
pub fn read_lambda(math: &Node) -> Result<(Vec<String>, Expr)> {
    let mut lambda = math.children.first().ok_or_else(|| bad("empty <math>"))?;
    if local(lambda) == "semantics" {
        lambda = lambda.children.first().ok_or_else(|| bad("empty semantics"))?;
    }
    if local(lambda) != "lambda" {
        return Err(bad("function definition without lambda"));
    }
    let mut args = Vec::new();
    let mut body = None;
    for child in &lambda.children {
        if local(child) == "bvar" {
            let ci = child.children.first().ok_or_else(|| bad("empty bvar"))?;
            args.push(ci.deep_text().trim().to_string());
        } else {
            body = Some(read_expr(child)?);
        }
    }
    Ok((args, body.ok_or_else(|| bad("lambda without body"))?))
}

fn local(node: &Node) -> &str {
    node.name.rsplit(':').next().unwrap_or(&node.name)
}

fn strip_prefix(node: &Node) -> &Node {
    node
}

fn read_expr(node: &Node) -> Result<Expr> {
    match local(node) {
        "ci" => Ok(Expr::Sym(node.deep_text().trim().to_string())),
        "cn" => read_cn(node),
        "true" => Ok(Expr::Num(1.0)),
        "false" => Ok(Expr::Num(0.0)),
        "pi" => Ok(Expr::Num(std::f64::consts::PI)),
        "exponentiale" => Ok(Expr::Num(std::f64::consts::E)),
        "infinity" => Ok(Expr::Num(f64::INFINITY)),
        "notanumber" => Ok(Expr::Num(f64::NAN)),
        "csymbol" => read_csymbol(node, &[]),
        "piecewise" => read_piecewise(node),
        "apply" => read_apply(node),
        "semantics" => read_expr(node.children.first().ok_or_else(|| bad("empty semantics"))?),
        other => Err(bad(format!("unsupported element <{other}>"))),
    }
}

fn read_cn(node: &Node) -> Result<Expr> {
    let parse = |s: &str| {
        let t = s.trim();
        match t {
            "INF" | "inf" | "Infinity" => Ok(f64::INFINITY),
            "-INF" | "-inf" | "-Infinity" => Ok(f64::NEG_INFINITY),
            "NaN" | "nan" => Ok(f64::NAN),
            _ => t.parse::<f64>().map_err(|_| bad(format!("bad number `{t}`"))),
        }
    };
    match node.attr("type") {
        Some("e-notation") => {
            let sep = node.children.first().ok_or_else(|| bad("e-notation without sep"))?;
            let text = format!("{}e{}", node.text.trim(), sep.tail.trim());
            Ok(Expr::Num(parse(&text)?))
        }
        Some("rational") => {
            let sep = node.children.first().ok_or_else(|| bad("rational without sep"))?;
            Ok(Expr::Num(parse(&node.text)? / parse(&sep.tail)?))
        }
        _ => Ok(Expr::Num(parse(&node.deep_text())?)),
    }
}

fn read_csymbol(node: &Node, args: &[Expr]) -> Result<Expr> {
    let url = node.attr("definitionURL").unwrap_or_default();
    if url.ends_with("/symbols/time") {
        Ok(Expr::Time)
    } else if url.ends_with("/symbols/avogadro") {
        Ok(Expr::Avogadro)
    } else if url.ends_with("/symbols/delay") {
        if args.len() != 2 {
            return Err(bad("delay requires two arguments"));
        }
        Ok(Expr::Delay(Box::new(args[0].clone()), Box::new(args[1].clone())))
    } else if url.ends_with("/symbols/rateOf") {
        match args {
            [Expr::Sym(name)] => Ok(Expr::RateOf(name.clone())),
            _ => Err(bad("rateOf requires a single identifier")),
        }
    } else {
        Err(bad(format!("unsupported csymbol `{url}`")))
    }
}

fn read_piecewise(node: &Node) -> Result<Expr> {
    let mut pieces = Vec::new();
    let mut otherwise = None;
    for child in &node.children {
        match local(child) {
            "piece" => {
                if child.children.len() != 2 {
                    return Err(bad("piece needs a value and a condition"));
                }
                pieces.push((read_expr(&child.children[0])?, read_expr(&child.children[1])?));
            }
            "otherwise" => {
                let inner = child.children.first().ok_or_else(|| bad("empty otherwise"))?;
                otherwise = Some(Box::new(read_expr(inner)?));
            }
            other => return Err(bad(format!("unexpected <{other}> in piecewise"))),
        }
    }
    Ok(Expr::Piecewise(pieces, otherwise))
}

fn op_named(name: &str) -> Option<Op> {
    Some(match name {
        "plus" => Op::Add,
        "minus" => Op::Sub,
        "times" => Op::Mul,
        "divide" => Op::Div,
        "power" => Op::Pow,
        "root" => Op::Root,
        "abs" => Op::Abs,
        "exp" => Op::Exp,
        "ln" => Op::Ln,
        "log" => Op::Log,
        "floor" => Op::Floor,
        "ceiling" => Op::Ceil,
        "factorial" => Op::Factorial,
        "eq" => Op::Eq,
        "neq" => Op::Neq,
        "gt" => Op::Gt,
        "lt" => Op::Lt,
        "geq" => Op::Geq,
        "leq" => Op::Leq,
        "and" => Op::And,
        "or" => Op::Or,
        "xor" => Op::Xor,
        "not" => Op::Not,
        "implies" => Op::Implies,
        "sin" => Op::Sin,
        "cos" => Op::Cos,
        "tan" => Op::Tan,
        "sec" => Op::Sec,
        "csc" => Op::Csc,
        "cot" => Op::Cot,
        "sinh" => Op::Sinh,
        "cosh" => Op::Cosh,
        "tanh" => Op::Tanh,
        "sech" => Op::Sech,
        "csch" => Op::Csch,
        "coth" => Op::Coth,
        "arcsin" => Op::Asin,
        "arccos" => Op::Acos,
        "arctan" => Op::Atan,
        "arcsec" => Op::Asec,
        "arccsc" => Op::Acsc,
        "arccot" => Op::Acot,
        "arcsinh" => Op::Asinh,
        "arccosh" => Op::Acosh,
        "arctanh" => Op::Atanh,
        "arcsech" => Op::Asech,
        "arccsch" => Op::Acsch,
        "arccoth" => Op::Acoth,
        "max" => Op::Max,
        "min" => Op::Min,
        "quotient" => Op::Quotient,
        "rem" => Op::Rem,
        _ => return None,
    })
}

fn read_apply(node: &Node) -> Result<Expr> {
    let head = node.children.first().ok_or_else(|| bad("empty apply"))?;
    let head_name = local(head);
    let mut qualifier: Option<Expr> = None;
    let mut args = Vec::new();
    for child in &node.children[1..] {
        match local(child) {
            "degree" | "logbase" => {
                let inner = child.children.first().ok_or_else(|| bad("empty qualifier"))?;
                qualifier = Some(read_expr(inner)?);
            }
            _ => args.push(read_expr(child)?),
        }
    }
    match head_name {
        "ci" => Ok(Expr::Call(head.deep_text().trim().to_string(), args)),
        "csymbol" => read_csymbol(head, &args),
        "root" => {
            let degree = qualifier.unwrap_or(Expr::Num(2.0));
            let mut full = vec![degree];
            full.extend(args);
            Ok(Expr::Apply(Op::Root, full))
        }
        "log" => {
            let base = qualifier.unwrap_or(Expr::Num(10.0));
            let mut full = vec![base];
            full.extend(args);
            Ok(Expr::Apply(Op::Log, full))
        }
        "minus" if args.len() == 1 => Ok(Expr::Apply(Op::Neg, args)),
        name => {
            let op = op_named(name).ok_or_else(|| bad(format!("unsupported operator <{name}>")))?;
            // n-ary plus and times with no or one operand have the identity value per MathML.
            match (op, args.len()) {
                (Op::Add, 0) => Ok(Expr::Num(0.0)),
                (Op::Mul, 0) => Ok(Expr::Num(1.0)),
                (Op::And, 0) => Ok(Expr::Num(1.0)),
                (Op::Or | Op::Xor, 0) => Ok(Expr::Num(0.0)),
                _ => Ok(Expr::Apply(op, args)),
            }
        }
    }
}

fn op_tag(op: Op) -> &'static str {
    match op {
        Op::Add => "plus",
        Op::Sub | Op::Neg => "minus",
        Op::Mul => "times",
        Op::Div => "divide",
        Op::Pow => "power",
        Op::Root => "root",
        Op::Abs => "abs",
        Op::Exp => "exp",
        Op::Ln => "ln",
        Op::Log => "log",
        Op::Floor => "floor",
        Op::Ceil => "ceiling",
        Op::Factorial => "factorial",
        Op::Eq => "eq",
        Op::Neq => "neq",
        Op::Gt => "gt",
        Op::Lt => "lt",
        Op::Geq => "geq",
        Op::Leq => "leq",
        Op::And => "and",
        Op::Or => "or",
        Op::Xor => "xor",
        Op::Not => "not",
        Op::Implies => "implies",
        Op::Sin => "sin",
        Op::Cos => "cos",
        Op::Tan => "tan",
        Op::Sec => "sec",
        Op::Csc => "csc",
        Op::Cot => "cot",
        Op::Sinh => "sinh",
        Op::Cosh => "cosh",
        Op::Tanh => "tanh",
        Op::Sech => "sech",
        Op::Csch => "csch",
        Op::Coth => "coth",
        Op::Asin => "arcsin",
        Op::Acos => "arccos",
        Op::Atan => "arctan",
        Op::Asec => "arcsec",
        Op::Acsc => "arccsc",
        Op::Acot => "arccot",
        Op::Asinh => "arcsinh",
        Op::Acosh => "arccosh",
        Op::Atanh => "arctanh",
        Op::Asech => "arcsech",
        Op::Acsch => "arccsch",
        Op::Acoth => "arccoth",
        Op::Max => "max",
        Op::Min => "min",
        Op::Quotient => "quotient",
        Op::Rem => "rem",
    }
}

/// Writes an expression as the content of a `<math>` element.
pub fn write_expr(e: &Expr, out: &mut String) {
    const CS: &str = "http://www.sbml.org/sbml/symbols";
    match e {
        Expr::Num(v) => {
            if v.is_nan() {
                out.push_str("<notanumber/>");
            } else if v.is_infinite() {
                if *v < 0.0 {
                    out.push_str("<apply><minus/><infinity/></apply>");
                } else {
                    out.push_str("<infinity/>");
                }
            } else {
                out.push_str(&format!("<cn>{v:?}</cn>"));
            }
        }
        Expr::Sym(s) => out.push_str(&format!("<ci>{}</ci>", super::write::escape(s))),
        Expr::Time => out.push_str(&format!("<csymbol encoding=\"text\" definitionURL=\"{CS}/time\">time</csymbol>")),
        Expr::Avogadro => {
            out.push_str(&format!("<csymbol encoding=\"text\" definitionURL=\"{CS}/avogadro\">avogadro</csymbol>"))
        }
        Expr::Apply(op, args) => {
            out.push_str(&format!("<apply><{}/>", op_tag(*op)));
            match op {
                Op::Root | Op::Log if args.len() == 2 => {
                    let q = if *op == Op::Root { "degree" } else { "logbase" };
                    out.push_str(&format!("<{q}>"));
                    write_expr(&args[0], out);
                    out.push_str(&format!("</{q}>"));
                    write_expr(&args[1], out);
                }
                _ => args.iter().for_each(|a| write_expr(a, out)),
            }
            out.push_str("</apply>");
        }
        Expr::Piecewise(pieces, otherwise) => {
            out.push_str("<piecewise>");
            for (value, cond) in pieces {
                out.push_str("<piece>");
                write_expr(value, out);
                write_expr(cond, out);
                out.push_str("</piece>");
            }
            if let Some(o) = otherwise {
                out.push_str("<otherwise>");
                write_expr(o, out);
                out.push_str("</otherwise>");
            }
            out.push_str("</piecewise>");
        }
        Expr::Call(name, args) => {
            out.push_str(&format!("<apply><ci>{}</ci>", super::write::escape(name)));
            args.iter().for_each(|a| write_expr(a, out));
            out.push_str("</apply>");
        }
        Expr::Delay(a, b) => {
            out.push_str(&format!("<apply><csymbol encoding=\"text\" definitionURL=\"{CS}/delay\">delay</csymbol>"));
            write_expr(a, out);
            write_expr(b, out);
            out.push_str("</apply>");
        }
        Expr::RateOf(s) => out.push_str(&format!(
            "<apply><csymbol encoding=\"text\" definitionURL=\"{CS}/rateOf\">rateOf</csymbol><ci>{}</ci></apply>",
            super::write::escape(s)
        )),
    }
}
