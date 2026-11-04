//! A small XML tree. SBML is deeply nested and read whole, so a tree keeps the reader simple.

use crate::error::{Error, Result};
use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};
use std::collections::HashMap;

/// An XML element. Names and attribute keys use canonical package prefixes (`fbc:`, `groups:`).
#[derive(Clone, Debug, Default)]
pub struct Node {
    /// Element name with canonical prefix.
    pub name: String,
    /// Attributes in document order.
    pub attrs: Vec<(String, String)>,
    /// Child elements.
    pub children: Vec<Node>,
    /// Text before the first child element.
    pub text: String,
    /// Text directly after this element's closing tag, inside the parent.
    pub tail: String,
}

impl Node {
    /// Value of an attribute.
    pub fn attr(&self, key: &str) -> Option<&str> {
        self.attrs.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
    }

    /// First child with the given name.
    pub fn child(&self, name: &str) -> Option<&Node> {
        self.children.iter().find(|c| c.name == name)
    }

    /// All children with the given name.
    pub fn children_named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Node> + 'a {
        self.children.iter().filter(move |c| c.name == name)
    }

    /// Text of this element and its descendants.
    pub fn deep_text(&self) -> String {
        let mut out = self.text.clone();
        for child in &self.children {
            out.push_str(&child.deep_text());
            out.push_str(&child.tail);
        }
        out
    }
}

fn canonical_package(uri: &str) -> Option<&'static str> {
    if uri.contains("/fbc/") {
        Some("fbc")
    } else if uri.contains("/groups/") {
        Some("groups")
    } else if uri.contains("/comp/") {
        Some("comp")
    } else if uri.contains("/layout/") {
        Some("layout")
    } else if uri.contains("/distrib/") {
        Some("distrib")
    } else if uri.contains("/render/") {
        Some("render")
    } else if uri.contains("/qual/") {
        Some("qual")
    } else {
        None
    }
}

struct Names {
    prefix_map: HashMap<String, String>,
}

impl Names {
    fn canonical(&self, raw: &str) -> String {
        match raw.split_once(':') {
            Some((prefix, local)) => match self.prefix_map.get(prefix) {
                Some(canon) => format!("{canon}:{local}"),
                None => format!("{prefix}:{local}"),
            },
            None => raw.to_string(),
        }
    }
}

fn open_node(start: &BytesStart<'_>, names: &mut Names, packages: &mut Vec<String>) -> Result<Node> {
    let bad = |e: &dyn std::fmt::Display| Error::Sbml(format!("xml: {e}"));
    let mut raw_attrs = Vec::new();
    for attr in start.attributes() {
        let attr = attr.map_err(|e| bad(&e))?;
        let key = std::str::from_utf8(attr.key.as_ref()).map_err(|e| bad(&e))?.to_string();
        let value = attr.unescape_value().map_err(|e| bad(&e))?.into_owned();
        if let Some(prefix) = key.strip_prefix("xmlns:")
            && let Some(canon) = canonical_package(&value)
        {
            names.prefix_map.insert(prefix.to_string(), canon.to_string());
            if !packages.iter().any(|p| p == canon) {
                packages.push(canon.to_string());
            }
        }
        raw_attrs.push((key, value));
    }
    let name = std::str::from_utf8(start.name().as_ref()).map_err(|e| bad(&e))?.to_string();
    Ok(Node {
        name: names.canonical(&name),
        attrs: raw_attrs.into_iter().map(|(k, v)| (names.canonical(&k), v)).collect(),
        ..Node::default()
    })
}

/// Parses an XML document into a tree and returns the root with the SBML packages it declares.
pub fn parse(xml: &str) -> Result<(Node, Vec<String>)> {
    let mut reader = Reader::from_str(xml);
    let mut names = Names { prefix_map: HashMap::new() };
    let mut packages = Vec::new();
    let mut stack: Vec<Node> = Vec::new();
    let mut root: Option<Node> = None;
    loop {
        let event = reader.read_event().map_err(|e| {
            Error::Sbml(format!("xml at byte {}: {e}", reader.buffer_position()))
        })?;
        match event {
            Event::Start(start) => stack.push(open_node(&start, &mut names, &mut packages)?),
            Event::Empty(start) => {
                let node = open_node(&start, &mut names, &mut packages)?;
                match stack.last_mut() {
                    Some(parent) => parent.children.push(node),
                    None => root = Some(node),
                }
            }
            Event::End(_) => {
                let node = stack.pop().ok_or_else(|| Error::Sbml("unbalanced xml".into()))?;
                match stack.last_mut() {
                    Some(parent) => parent.children.push(node),
                    None => root = Some(node),
                }
            }
            Event::Text(text) => {
                if let Some(parent) = stack.last_mut() {
                    let raw = text.decode().map_err(|e| Error::Sbml(format!("xml: {e}")))?;
                    let value = quick_xml::escape::unescape(&raw)
                        .map_err(|e| Error::Sbml(format!("xml: {e}")))?;
                    append_text(parent, &value);
                }
            }
            Event::GeneralRef(entity) => {
                if let Some(parent) = stack.last_mut() {
                    let raw = entity.decode().map_err(|e| Error::Sbml(format!("xml: {e}")))?;
                    let value = quick_xml::escape::unescape(&format!("&{raw};"))
                        .map_err(|e| Error::Sbml(format!("xml: {e}")))?
                        .into_owned();
                    append_text(parent, &value);
                }
            }
            Event::CData(text) => {
                if let Some(parent) = stack.last_mut() {
                    append_text(parent, &String::from_utf8_lossy(&text));
                }
            }
            Event::Eof => break,
            _ => {}
        }
    }
    let root = root.ok_or_else(|| Error::Sbml("empty document".into()))?;
    Ok((root, packages))
}

fn append_text(parent: &mut Node, value: &str) {
    match parent.children.last_mut() {
        Some(last) => last.tail.push_str(value),
        None => parent.text.push_str(value),
    }
}
