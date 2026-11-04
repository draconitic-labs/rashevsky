//! Gene-protein-reaction rules: Boolean formulas over genes.

use crate::error::{Error, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A gene-protein-reaction rule. `And` joins subunits of a complex, `Or` joins isozymes.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Gpr {
    /// No genes are associated; the reaction cannot be removed by a gene deletion.
    #[default]
    None,
    /// A single gene, by index into the model's gene list.
    Gene(usize),
    /// All operands are required.
    And(Vec<Gpr>),
    /// Any operand suffices.
    Or(Vec<Gpr>),
}

impl Gpr {
    /// True if the reaction can still be catalysed when the genes in `deleted` are absent.
    pub fn active(&self, deleted: &[bool]) -> bool {
        match self {
            Gpr::None => true,
            Gpr::Gene(g) => !deleted[*g],
            Gpr::And(items) => items.iter().all(|g| g.active(deleted)),
            Gpr::Or(items) => items.iter().any(|g| g.active(deleted)),
        }
    }

    /// Gene indices that appear in the rule.
    pub fn genes(&self, out: &mut Vec<usize>) {
        match self {
            Gpr::None => {}
            Gpr::Gene(g) => out.push(*g),
            Gpr::And(items) | Gpr::Or(items) => items.iter().for_each(|g| g.genes(out)),
        }
    }

    /// Text form, with `name` mapping gene indices to identifiers.
    pub fn to_rule(&self, name: &impl Fn(usize) -> String) -> String {
        match self {
            Gpr::None => String::new(),
            Gpr::Gene(g) => name(*g),
            Gpr::And(items) => join(items, " and ", name),
            Gpr::Or(items) => join(items, " or ", name),
        }
    }

    /// Parses a rule such as `b0001 and (b0002 or b0003)`, adding unknown genes through `intern`.
    pub fn parse(text: &str, intern: &mut impl FnMut(&str) -> usize) -> Result<Self> {
        let tokens = tokenize(text);
        if tokens.is_empty() {
            return Ok(Gpr::None);
        }
        let mut parser = Parser { tokens: &tokens, pos: 0, intern };
        let gpr = parser.or_expr()?;
        if parser.pos != tokens.len() {
            return Err(Error::Gpr(format!("trailing input in `{text}`")));
        }
        Ok(gpr)
    }
}

fn join(items: &[Gpr], sep: &str, name: &impl Fn(usize) -> String) -> String {
    let mut parts = Vec::with_capacity(items.len());
    for item in items {
        let text = item.to_rule(name);
        parts.push(match item {
            Gpr::And(_) | Gpr::Or(_) => format!("({text})"),
            _ => text,
        });
    }
    parts.join(sep)
}

#[derive(Clone, Debug, PartialEq)]
enum Token {
    And,
    Or,
    Open,
    Close,
    Id(String),
}

fn tokenize(text: &str) -> Vec<Token> {
    let spaced = text.replace('(', " ( ").replace(')', " ) ").replace('&', " & ").replace('|', " | ");
    let mut tokens: Vec<Token> = spaced
        .split_whitespace()
        .map(|w| match w {
            "(" => Token::Open,
            ")" => Token::Close,
            _ if w.eq_ignore_ascii_case("and") || w == "&" => Token::And,
            _ if w.eq_ignore_ascii_case("or") || w == "|" => Token::Or,
            _ => Token::Id(w.to_string()),
        })
        .collect();
    // `&&` and `||` arrive as doubled operators.
    tokens.dedup_by(|a, b| matches!((a, b), (Token::And, Token::And) | (Token::Or, Token::Or)));
    tokens
}

struct Parser<'a, F: FnMut(&str) -> usize> {
    tokens: &'a [Token],
    pos: usize,
    intern: &'a mut F,
}

impl<F: FnMut(&str) -> usize> Parser<'_, F> {
    fn or_expr(&mut self) -> Result<Gpr> {
        let mut items = vec![self.and_expr()?];
        while self.tokens.get(self.pos) == Some(&Token::Or) {
            self.pos += 1;
            items.push(self.and_expr()?);
        }
        Ok(if items.len() == 1 { items.remove(0) } else { Gpr::Or(items) })
    }

    fn and_expr(&mut self) -> Result<Gpr> {
        let mut items = vec![self.atom()?];
        while self.tokens.get(self.pos) == Some(&Token::And) {
            self.pos += 1;
            items.push(self.atom()?);
        }
        Ok(if items.len() == 1 { items.remove(0) } else { Gpr::And(items) })
    }

    fn atom(&mut self) -> Result<Gpr> {
        match self.tokens.get(self.pos) {
            Some(Token::Open) => {
                self.pos += 1;
                let inner = self.or_expr()?;
                if self.tokens.get(self.pos) != Some(&Token::Close) {
                    return Err(Error::Gpr("missing closing parenthesis".into()));
                }
                self.pos += 1;
                Ok(inner)
            }
            Some(Token::Id(id)) => {
                self.pos += 1;
                Ok(Gpr::Gene((self.intern)(id)))
            }
            other => Err(Error::Gpr(format!("unexpected token {other:?}"))),
        }
    }
}

/// Index of gene identifiers used while parsing rules.
#[derive(Default, Debug)]
pub struct GeneIndex {
    /// Identifier to position.
    pub by_id: HashMap<String, usize>,
    /// Identifiers in order of first use.
    pub ids: Vec<String>,
}

impl GeneIndex {
    /// Returns the index of `id`, adding it if new.
    pub fn intern(&mut self, id: &str) -> usize {
        if let Some(&i) = self.by_id.get(id) {
            return i;
        }
        let i = self.ids.len();
        self.ids.push(id.to_string());
        self.by_id.insert(id.to_string(), i);
        i
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(text: &str) -> (Gpr, GeneIndex) {
        let mut idx = GeneIndex::default();
        let gpr = Gpr::parse(text, &mut |s| idx.intern(s)).unwrap();
        (gpr, idx)
    }

    #[test]
    fn and_binds_tighter_than_or() {
        let (gpr, idx) = parse("a or b and c");
        let del = |names: &[&str]| {
            let mut v = vec![false; idx.ids.len()];
            for n in names {
                v[idx.by_id[*n]] = true;
            }
            v
        };
        assert!(gpr.active(&del(&["b"])));
        assert!(gpr.active(&del(&["a"])));
        assert!(!gpr.active(&del(&["a", "b"])));
        assert!(gpr.active(&del(&["b", "c"])));
    }

    #[test]
    fn empty_rule_cannot_be_deleted() {
        let (gpr, _) = parse("  ");
        assert_eq!(gpr, Gpr::None);
        assert!(gpr.active(&[]));
    }

    #[test]
    fn text_form_round_trips() {
        let (gpr, idx) = parse("(a and b) or (c and (d or e))");
        let text = gpr.to_rule(&|g| idx.ids[g].clone());
        let mut idx2 = GeneIndex::default();
        let gpr2 = Gpr::parse(&text, &mut |s| idx2.intern(s)).unwrap();
        assert_eq!(gpr, gpr2);
    }

    #[test]
    fn rejects_unbalanced_parentheses() {
        let mut idx = GeneIndex::default();
        assert!(Gpr::parse("(a or b", &mut |s| idx.intern(s)).is_err());
        assert!(Gpr::parse("a or", &mut |s| idx.intern(s)).is_err());
    }
}
