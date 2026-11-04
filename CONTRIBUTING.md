# Contributing

Run these before sending a change:

```
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --release
maturin develop --release && pytest python/tests
```

A change to a numerical method needs a test that fails without it, and an update to `docs/MODEL.md`.
Run `python validation/run_all.py` when a change could move a validation number. Never change a
tolerance to make a benchmark pass: fix the engine, or explain the benchmark's limit in
`docs/VALIDATION.md` next to the unchanged number.

The desktop interface is tested in its own Tauri window, with `desktop/scripts/qa.py`, and not in a
browser: the webview it ships in behaves differently.
