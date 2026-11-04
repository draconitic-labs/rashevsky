# Labs demo

`e_coli_core` flux balance in the browser, on the pure Rust simplex. The full engine is the product.

Build with `./build.sh`, which needs the `wasm32-unknown-unknown` target and `wasm-bindgen-cli` at the
version pinned in `crates/rashevsky-wasm/Cargo.toml`. Serve `web/dist` as static files, with a content
security policy that allows `script-src 'self' 'wasm-unsafe-eval'`.

Test in the webview it targets, not a browser:

```
(cd web/preview-shell && cargo build --release)
```

then run `rashevsky-web-preview` from the target folder, or drive it with `desktop/scripts/drive.py`.

The demo bundles `e_coli_core` from BiGG Models, which is licensed for educational, research and
non-profit use. See `models/e_coli_core/NOTICE` before placing it on a commercial site.
