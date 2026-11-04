//! Case files: a study described in JSON and run with `rashevsky run`.

use rashevsky::Result;
use std::path::{Path, PathBuf};

/// Runs the study in `path`.
pub fn run(path: &Path, _out: &Option<PathBuf>) -> Result<()> {
    Err(rashevsky::Error::Invalid(format!("case file {} is not supported yet", path.display())))
}
