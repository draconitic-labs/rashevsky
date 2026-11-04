//! HiGHS backend over its C API.
//!
//! Each session owns one `Highs` instance running on one thread. Parallel analyses open one session
//! per worker, so the thread pool of the engine and not that of HiGHS decides how work is shared.

use super::{Backend, Problem, Sense, Session, Solution, Status};
use crate::error::{Error, Result};
use highs_sys as ffi;
use std::ffi::CString;
use std::os::raw::{c_int, c_void};

/// HiGHS solver settings.
#[derive(Clone, Debug)]
pub struct Highs {
    /// Primal and dual feasibility tolerance.
    pub feasibility_tolerance: f64,
    /// Relative MIP gap at which a mixed-integer solve stops.
    pub mip_rel_gap: f64,
    /// Integrality tolerance for mixed-integer solves.
    pub mip_feasibility_tolerance: f64,
    /// Wall-clock limit in seconds per solve. Infinite by default.
    pub time_limit: f64,
}

impl Default for Highs {
    fn default() -> Self {
        Self {
            feasibility_tolerance: 1e-9,
            mip_rel_gap: 1e-6,
            mip_feasibility_tolerance: 1e-9,
            time_limit: f64::INFINITY,
        }
    }
}

impl Backend for Highs {
    fn name(&self) -> &'static str {
        "highs"
    }

    fn open(&self, problem: &Problem) -> Result<Box<dyn Session>> {
        Ok(Box::new(HighsSession::new(problem, self)?))
    }
}

struct HighsSession {
    ptr: *mut c_void,
    n_cols: usize,
    n_rows: usize,
    is_mip: bool,
    is_qp: bool,
    hessian: Vec<(usize, usize, f64)>,
    obj: Vec<f64>,
    sense: Sense,
}

// A HiGHS instance is used by one thread at a time; the pointer moves with the session.
unsafe impl Send for HighsSession {}

impl Drop for HighsSession {
    fn drop(&mut self) {
        // SAFETY: `ptr` came from `Highs_create` and is destroyed exactly once.
        unsafe { ffi::Highs_destroy(self.ptr) }
    }
}

const HESSIAN_FORMAT_TRIANGULAR: c_int = 1;

fn check(code: c_int, what: &str) -> Result<()> {
    if code == ffi::STATUS_ERROR {
        return Err(Error::Solver(format!("HiGHS call failed: {what}")));
    }
    Ok(())
}

fn inf_clamp(v: f64) -> f64 {
    v.clamp(f64::NEG_INFINITY, f64::INFINITY)
}

impl HighsSession {
    fn new(p: &Problem, opts: &Highs) -> Result<Self> {
        // SAFETY: `Highs_create` has no preconditions.
        let ptr = unsafe { ffi::Highs_create() };
        if ptr.is_null() {
            return Err(Error::Solver("Highs_create returned null".into()));
        }
        let sense = p.sense.unwrap_or(Sense::Minimize);
        let mut session = Self {
            ptr,
            n_cols: p.n_cols(),
            n_rows: p.n_rows(),
            is_mip: p.is_mip(),
            is_qp: p.is_qp(),
            hessian: p.hessian.clone(),
            obj: p.obj.clone(),
            sense,
        };
        session.set_bool("output_flag", false)?;
        session.set_int("threads", 1)?;
        session.set_double("primal_feasibility_tolerance", opts.feasibility_tolerance)?;
        session.set_double("dual_feasibility_tolerance", opts.feasibility_tolerance)?;
        session.set_double("mip_rel_gap", opts.mip_rel_gap)?;
        session.set_double("mip_abs_gap", 1e-9)?;
        session.set_double("mip_feasibility_tolerance", opts.mip_feasibility_tolerance)?;
        if opts.time_limit.is_finite() {
            session.set_double("time_limit", opts.time_limit)?;
        }
        session.pass(p)?;
        Ok(session)
    }

    fn set_bool(&mut self, name: &str, value: bool) -> Result<()> {
        let key = CString::new(name).expect("option name");
        // SAFETY: valid instance and NUL-terminated option name.
        let code = unsafe { ffi::Highs_setBoolOptionValue(self.ptr, key.as_ptr(), value as c_int) };
        check(code, name)
    }

    fn set_int(&mut self, name: &str, value: c_int) -> Result<()> {
        let key = CString::new(name).expect("option name");
        // SAFETY: valid instance and NUL-terminated option name.
        let code = unsafe { ffi::Highs_setIntOptionValue(self.ptr, key.as_ptr(), value) };
        check(code, name)
    }

    fn set_double(&mut self, name: &str, value: f64) -> Result<()> {
        let key = CString::new(name).expect("option name");
        // SAFETY: valid instance and NUL-terminated option name.
        let code = unsafe { ffi::Highs_setDoubleOptionValue(self.ptr, key.as_ptr(), value) };
        check(code, name)
    }

    fn pass(&mut self, p: &Problem) -> Result<()> {
        let (start, index, value) = p.to_csc();
        let start: Vec<c_int> = start.iter().map(|&v| v as c_int).collect();
        let index: Vec<c_int> = index.iter().map(|&v| v as c_int).collect();
        let col_lb: Vec<f64> = p.col_lb.iter().map(|&v| inf_clamp(v)).collect();
        let col_ub: Vec<f64> = p.col_ub.iter().map(|&v| inf_clamp(v)).collect();
        let row_lb: Vec<f64> = p.row_lb.iter().map(|&v| inf_clamp(v)).collect();
        let row_ub: Vec<f64> = p.row_ub.iter().map(|&v| inf_clamp(v)).collect();
        let sense = match self.sense {
            Sense::Minimize => ffi::OBJECTIVE_SENSE_MINIMIZE,
            Sense::Maximize => ffi::OBJECTIVE_SENSE_MAXIMIZE,
        };
        let integrality: Vec<c_int> = if p.is_mip() {
            p.integer.iter().map(|&b| if b { 1 } else { 0 }).collect()
        } else {
            Vec::new()
        };
        // SAFETY: every pointer refers to a live vector of the length HiGHS is told; an empty
        // integrality array is passed as null, which HiGHS reads as all continuous.
        let code = unsafe {
            ffi::Highs_passModel(
                self.ptr,
                p.n_cols() as c_int,
                p.n_rows() as c_int,
                index.len() as c_int,
                0,
                ffi::MATRIX_FORMAT_COLUMN_WISE,
                HESSIAN_FORMAT_TRIANGULAR,
                sense,
                0.0,
                p.obj.as_ptr(),
                col_lb.as_ptr(),
                col_ub.as_ptr(),
                row_lb.as_ptr(),
                row_ub.as_ptr(),
                start.as_ptr(),
                index.as_ptr(),
                value.as_ptr(),
                std::ptr::null(),
                std::ptr::null(),
                std::ptr::null(),
                if integrality.is_empty() {
                    std::ptr::null()
                } else {
                    integrality.as_ptr()
                },
            )
        };
        check(code, "passModel")?;
        if p.is_qp() {
            self.pass_hessian(p)?;
        }
        Ok(())
    }

    fn pass_hessian(&mut self, p: &Problem) -> Result<()> {
        // HiGHS wants the lower triangle column-wise.
        let n = p.n_cols();
        let mut entries = p.hessian.clone();
        entries.sort_by_key(|&(i, j, _)| (j, i));
        let mut start = vec![0 as c_int; n + 1];
        let mut index = Vec::with_capacity(entries.len());
        let mut value = Vec::with_capacity(entries.len());
        for &(i, j, q) in &entries {
            if i < j {
                return Err(Error::Solver("hessian entries must have row >= col".into()));
            }
            start[j + 1] += 1;
            index.push(i as c_int);
            value.push(q);
        }
        for j in 0..n {
            start[j + 1] += start[j];
        }
        // SAFETY: pointers refer to live vectors; lengths match the arguments.
        let code = unsafe {
            ffi::Highs_passHessian(
                self.ptr,
                n as c_int,
                value.len() as c_int,
                HESSIAN_FORMAT_TRIANGULAR,
                start.as_ptr(),
                index.as_ptr(),
                value.as_ptr(),
            )
        };
        check(code, "passHessian")
    }
}

impl Session for HighsSession {
    fn set_col_bounds(&mut self, col: usize, lb: f64, ub: f64) -> Result<()> {
        // SAFETY: valid instance; the column index is checked by HiGHS.
        let code = unsafe { ffi::Highs_changeColBounds(self.ptr, col as c_int, lb, ub) };
        check(code, "changeColBounds")
    }

    fn set_row_bounds(&mut self, row: usize, lb: f64, ub: f64) -> Result<()> {
        // SAFETY: valid instance; the row index is checked by HiGHS.
        let code = unsafe { ffi::Highs_changeRowBounds(self.ptr, row as c_int, lb, ub) };
        check(code, "changeRowBounds")
    }

    fn set_objective(&mut self, sense: Sense, obj: &[f64]) -> Result<()> {
        if obj.len() != self.n_cols {
            return Err(Error::Solver("objective length differs from column count".into()));
        }
        if sense != self.sense {
            let s = match sense {
                Sense::Minimize => ffi::OBJECTIVE_SENSE_MINIMIZE,
                Sense::Maximize => ffi::OBJECTIVE_SENSE_MAXIMIZE,
            };
            // SAFETY: valid instance.
            check(unsafe { ffi::Highs_changeObjectiveSense(self.ptr, s) }, "sense")?;
            self.sense = sense;
        }
        for (j, (&new, old)) in obj.iter().zip(self.obj.iter_mut()).enumerate() {
            if new != *old {
                // SAFETY: valid instance; `j` is within the column count checked above.
                check(unsafe { ffi::Highs_changeColCost(self.ptr, j as c_int, new) }, "cost")?;
                *old = new;
            }
        }
        Ok(())
    }

    fn add_row(&mut self, lb: f64, ub: f64, coefs: &[(usize, f64)]) -> Result<usize> {
        let index: Vec<c_int> = coefs.iter().map(|&(j, _)| j as c_int).collect();
        let value: Vec<f64> = coefs.iter().map(|&(_, v)| v).collect();
        // SAFETY: pointers refer to live vectors of the stated length.
        let code = unsafe {
            ffi::Highs_addRow(self.ptr, lb, ub, index.len() as c_int, index.as_ptr(), value.as_ptr())
        };
        check(code, "addRow")?;
        self.n_rows += 1;
        Ok(self.n_rows - 1)
    }

    fn solve(&mut self) -> Result<Solution> {
        // SAFETY: valid instance.
        let run = unsafe { ffi::Highs_run(self.ptr) };
        check(run, "run")?;
        // SAFETY: valid instance.
        let model_status = unsafe { ffi::Highs_getModelStatus(self.ptr) };
        let status = match model_status {
            ffi::MODEL_STATUS_OPTIMAL | ffi::MODEL_STATUS_MODEL_EMPTY => Status::Optimal,
            ffi::MODEL_STATUS_INFEASIBLE => Status::Infeasible,
            ffi::MODEL_STATUS_UNBOUNDED => Status::Unbounded,
            ffi::MODEL_STATUS_UNBOUNDED_OR_INFEASIBLE => self.classify_unbounded_or_infeasible()?,
            _ => Status::Unknown,
        };
        let mut x = vec![0.0; self.n_cols];
        let mut col_dual = vec![0.0; self.n_cols];
        let mut row_activity = vec![0.0; self.n_rows];
        let mut row_dual = vec![0.0; self.n_rows];
        // SAFETY: each buffer has the size HiGHS writes (columns, then rows).
        let code = unsafe {
            ffi::Highs_getSolution(
                self.ptr,
                x.as_mut_ptr(),
                col_dual.as_mut_ptr(),
                row_activity.as_mut_ptr(),
                row_dual.as_mut_ptr(),
            )
        };
        check(code, "getSolution")?;
        // SAFETY: valid instance.
        let objective = unsafe { ffi::Highs_getObjectiveValue(self.ptr) };
        let (row_dual, col_dual) = if self.is_mip {
            (Vec::new(), Vec::new())
        } else {
            (row_dual, col_dual)
        };
        let _ = (&self.hessian, self.is_qp);
        Ok(Solution {
            status,
            objective: if status == Status::Optimal { objective } else { f64::NAN },
            x,
            row_activity,
            row_dual,
            col_dual,
        })
    }
}

impl HighsSession {
    /// Presolve can stop at "unbounded or infeasible". Re-solving with presolve off decides which.
    fn classify_unbounded_or_infeasible(&mut self) -> Result<Status> {
        let key = CString::new("presolve").expect("option");
        let off = CString::new("off").expect("value");
        let on = CString::new("choose").expect("value");
        // SAFETY: valid instance and NUL-terminated strings.
        check(unsafe { ffi::Highs_setStringOptionValue(self.ptr, key.as_ptr(), off.as_ptr()) }, "presolve")?;
        // SAFETY: valid instance.
        check(unsafe { ffi::Highs_run(self.ptr) }, "run")?;
        // SAFETY: valid instance.
        let status = unsafe { ffi::Highs_getModelStatus(self.ptr) };
        // SAFETY: valid instance and NUL-terminated strings.
        check(unsafe { ffi::Highs_setStringOptionValue(self.ptr, key.as_ptr(), on.as_ptr()) }, "presolve")?;
        Ok(match status {
            ffi::MODEL_STATUS_INFEASIBLE => Status::Infeasible,
            ffi::MODEL_STATUS_UNBOUNDED => Status::Unbounded,
            ffi::MODEL_STATUS_OPTIMAL => Status::Optimal,
            _ => Status::Unknown,
        })
    }
}
