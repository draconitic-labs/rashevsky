## What changes

## How it was checked

- [ ] `cargo fmt --all --check` and `cargo clippy --workspace --all-targets -- -D warnings`
- [ ] `cargo test --workspace --release`
- [ ] A change to a numerical method has a test that fails without it, and an update to `docs/MODEL.md`
- [ ] `python validation/run_all.py`, if a validation number could move
