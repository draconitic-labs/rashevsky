//! Rashevsky: genome-scale constraint-based and detailed kinetic metabolic modelling in one library.
//!
//! The equations, numerical methods and references are in `docs/MODEL.md`; the validation numbers
//! and the commands that reproduce them are in `docs/VALIDATION.md`.

#![warn(missing_docs)]

pub mod cb;
pub mod error;
pub mod gpr;
pub mod io;
pub mod kinetic;
pub mod linalg;
pub mod lp;
pub mod model;
pub mod parallel;
pub mod sbml;
pub mod suite;
pub mod tasks;

pub use error::{Error, Result};
pub use nalgebra;
pub use model::Model;
