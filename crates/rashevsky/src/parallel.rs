//! Data parallelism over the Rayon pool, with a sequential fallback for WebAssembly builds.

/// Work unit size: a worker keeps one solver session per chunk, and a session's first solve is much
/// the slowest, so chunks stay large while still leaving room for balancing between workers.
pub fn chunk_size(n: usize) -> usize {
    #[cfg(feature = "parallel")]
    let workers = rayon::current_num_threads().max(1);
    #[cfg(not(feature = "parallel"))]
    let workers = 1usize;
    (n / (workers * 2)).max(1)
}

/// Applies `f` to chunks of `items`, keeping the order of chunks.
pub fn map_chunks<T, R, F>(items: &[T], chunk: usize, f: F) -> Vec<R>
where
    T: Sync,
    R: Send,
    F: Fn(&[T]) -> R + Sync + Send,
{
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        items.par_chunks(chunk.max(1)).map(f).collect()
    }
    #[cfg(not(feature = "parallel"))]
    {
        items.chunks(chunk.max(1)).map(f).collect()
    }
}
