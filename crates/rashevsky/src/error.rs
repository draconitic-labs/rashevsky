//! Error type shared by every module of the engine.

use thiserror::Error;

/// Errors returned by the engine.
#[derive(Debug, Error)]
pub enum Error {
    /// A file could not be read or written.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    /// The SBML document is malformed or uses a construct the reader does not accept.
    #[error("sbml: {0}")]
    Sbml(String),
    /// A gene-protein-reaction rule could not be parsed.
    #[error("gpr: {0}")]
    Gpr(String),
    /// The optimizer failed or returned a status other than the one required.
    #[error("solver: {0}")]
    Solver(String),
    /// A kinetic model could not be compiled or simulated.
    #[error("kinetic: {0}")]
    Kinetic(String),
    /// The numerical integrator failed.
    #[error("integrator: {0}")]
    Integrator(String),
    /// An identifier does not name an object in the model.
    #[error("unknown {kind} `{id}`")]
    Unknown {
        /// Kind of object that was looked up.
        kind: &'static str,
        /// Identifier that was not found.
        id: String,
    },
    /// The arguments of an analysis are inconsistent.
    #[error("invalid input: {0}")]
    Invalid(String),
    /// A JSON document could not be read or written.
    #[error("json: {0}")]
    Json(#[from] serde_json::Error),
}

/// Result alias used across the crate.
pub type Result<T> = std::result::Result<T, Error>;
