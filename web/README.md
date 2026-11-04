# Labs demo

iCH360 flux balance in the browser, on the pure Rust simplex. The full engine is the product.

Build with `./build.sh`, which needs the `wasm32-unknown-unknown` target and `wasm-bindgen-cli` at the
version pinned in `crates/rashevsky-wasm/Cargo.toml`. Serve `web/dist` as static files, with a content
security policy that allows `script-src 'self' 'wasm-unsafe-eval'`.

Test in the webview it targets, not a browser:

```
(cd web/preview-shell && cargo build --release)
```

then run `rashevsky-web-preview` from the target folder, or drive it with `desktop/scripts/drive.py`.

The demo bundles iCH360. Its notice has to travel with the model: `build.sh` copies it to
`dist/NOTICE-iCH360.txt` and the licenses of linked software to `dist/THIRD-PARTY-NOTICES.txt`, and the
page links to both. Keep them when installing. Before putting the demo on a commercial site, read the
question about the model in `docs/DATA.md`.
