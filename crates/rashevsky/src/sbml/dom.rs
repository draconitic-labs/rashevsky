//! A small XML tree. SBML is deeply nested and read whole, so a tree keeps the reader simple.

use crate::error::{Error, Result};
use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};

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

/// Namespace bindings in force, innermost last. A binding lasts as long as its element.
struct Names {
    /// `(element depth, prefix, canonical name)`, with an empty prefix for the default namespace.
    bindings: Vec<(usize, String, String)>,
    depth: usize,
}

impl Names {
    fn lookup(&self, prefix: &str) -> Option<&str> {
        self.bindings.iter().rev().find(|(_, p, _)| p == prefix).map(|(_, _, canon)| canon.as_str())
    }

    fn close(&mut self) {
        let depth = self.depth;
        self.bindings.retain(|(d, _, _)| *d < depth);
        self.depth -= 1;
    }

    fn canonical(&self, raw: &str, is_element: bool) -> String {
        match raw.split_once(':') {
            Some((prefix, local)) => match self.lookup(prefix) {
                // MathML and SBML core have no prefix once their namespaces are known.
                Some("") => local.to_string(),
                Some(canon) => format!("{canon}:{local}"),
                None => format!("{prefix}:{local}"),
            },
            // A default namespace applies to elements only, never to attributes.
            None if is_element => match self.lookup("") {
                Some(canon) if !canon.is_empty() => format!("{canon}:{raw}"),
                _ => raw.to_string(),
            },
            None => raw.to_string(),
        }
    }
}

/// The name a namespace goes by in the tree: a package name or `rdf`, nothing for MathML and SBML core, and
/// for a namespace this reader does not interpret, the prefix as written.
fn canonical_namespace(uri: &str, prefix: &str) -> String {
    if let Some(package) = canonical_package(uri) {
        return package.to_string();
    }
    if uri.contains("22-rdf-syntax-ns") {
        return "rdf".to_string();
    }
    let core = uri.contains("sbml.org/sbml/level") && (uri.ends_with("/core") || !uri.contains("/level3/"));
    if core || uri.contains("/Math/MathML") { String::new() } else { prefix.to_string() }
}

fn open_node(start: &BytesStart<'_>, names: &mut Names, packages: &mut Vec<String>) -> Result<Node> {
    let bad = |e: &dyn std::fmt::Display| Error::Sbml(format!("xml: {e}"));
    names.depth += 1;
    let mut raw_attrs = Vec::new();
    for attr in start.attributes() {
        let attr = attr.map_err(|e| bad(&e))?;
        let key = std::str::from_utf8(attr.key.as_ref()).map_err(|e| bad(&e))?.to_string();
        let value = attr.unescape_value().map_err(|e| bad(&e))?.into_owned();
        let prefix = if key == "xmlns" { Some("") } else { key.strip_prefix("xmlns:") };
        if let Some(prefix) = prefix {
            names.bindings.push((names.depth, prefix.to_string(), canonical_namespace(&value, prefix)));
            if let Some(canon) = canonical_package(&value)
                && !packages.iter().any(|p| p == canon)
            {
                packages.push(canon.to_string());
            }
        }
        raw_attrs.push((key, value));
    }
    let name = std::str::from_utf8(start.name().as_ref()).map_err(|e| bad(&e))?.to_string();
    Ok(Node {
        name: names.canonical(&name, true),
        attrs: raw_attrs.into_iter().map(|(k, v)| (names.canonical(&k, false), v)).collect(),
        ..Node::default()
    })
}

/// Parses an XML document into a tree and returns the root with the SBML packages it declares.
pub fn parse(xml: &str) -> Result<(Node, Vec<String>)> {
    let mut reader = Reader::from_str(xml);
    let mut names = Names { bindings: Vec::new(), depth: 0 };
    let mut packages = Vec::new();
    let mut stack: Vec<Node> = Vec::new();
    let mut root: Option<Node> = None;
    loop {
        let event =
            reader.read_event().map_err(|e| Error::Sbml(format!("xml at byte {}: {e}", reader.buffer_position())))?;
        match event {
            Event::Start(start) => stack.push(open_node(&start, &mut names, &mut packages)?),
            Event::Empty(start) => {
                let node = open_node(&start, &mut names, &mut packages)?;
                names.close();
                match stack.last_mut() {
                    Some(parent) => parent.children.push(node),
                    None => root = Some(node),
                }
            }
            Event::End(_) => {
                let node = stack.pop().ok_or_else(|| Error::Sbml("unbalanced xml".into()))?;
                names.close();
                match stack.last_mut() {
                    Some(parent) => parent.children.push(node),
                    None => root = Some(node),
                }
            }
            Event::Text(text) => {
                if let Some(parent) = stack.last_mut() {
                    let raw = text.decode().map_err(|e| Error::Sbml(format!("xml: {e}")))?;
                    let value = quick_xml::escape::unescape(&raw).map_err(|e| Error::Sbml(format!("xml: {e}")))?;
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
