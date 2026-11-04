//! Long analyses run on background threads; the interface polls them and can cancel them.

use serde::Serialize;
use serde_json::Value;
use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

/// A running or finished analysis.
pub struct Job {
    pub cancelled: AtomicBool,
    pub done: AtomicUsize,
    pub total: AtomicUsize,
    pub result: Mutex<Option<Result<Value, String>>>,
    pub started: Instant,
}

impl Job {
    pub fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }
}

/// All jobs of this session.
#[derive(Default)]
pub struct Jobs {
    next: AtomicU64,
    items: Mutex<HashMap<u64, Arc<Job>>>,
}

/// What a poll returns.
#[derive(Serialize)]
pub struct JobStatus {
    pub done: usize,
    pub total: usize,
    pub seconds: f64,
    pub finished: bool,
    pub result: Option<Value>,
    pub error: Option<String>,
}

impl Jobs {
    pub fn spawn(&self, total: usize, work: impl FnOnce(&Job) -> Result<Value, String> + Send + 'static) -> u64 {
        let id = self.next.fetch_add(1, Ordering::Relaxed);
        let job = Arc::new(Job {
            cancelled: AtomicBool::new(false),
            done: AtomicUsize::new(0),
            total: AtomicUsize::new(total),
            result: Mutex::new(None),
            started: Instant::now(),
        });
        self.items.lock().expect("jobs").insert(id, job.clone());
        std::thread::spawn(move || {
            let outcome = work(&job);
            *job.result.lock().expect("result") = Some(outcome);
        });
        id
    }

    pub fn poll(&self, id: u64) -> Option<JobStatus> {
        let job = self.items.lock().expect("jobs").get(&id)?.clone();
        let taken = job.result.lock().expect("result").take();
        let finished = taken.is_some();
        if finished {
            self.items.lock().expect("jobs").remove(&id);
        }
        let (result, error) = match taken {
            Some(Ok(v)) => (Some(v), None),
            Some(Err(e)) => (None, Some(e)),
            None => (None, None),
        };
        Some(JobStatus {
            done: job.done.load(Ordering::Relaxed),
            total: job.total.load(Ordering::Relaxed),
            seconds: job.started.elapsed().as_secs_f64(),
            finished,
            result,
            error,
        })
    }

    pub fn cancel(&self, id: u64) {
        if let Some(job) = self.items.lock().expect("jobs").get(&id) {
            job.cancelled.store(true, Ordering::Relaxed);
        }
    }
}
