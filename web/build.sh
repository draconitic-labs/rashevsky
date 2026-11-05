#!/usr/bin/env bash
# Builds the Labs demo into web/dist: WebAssembly, its JavaScript glue and the page.
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
root="$here/.."
target="${CARGO_TARGET_DIR:-$root/target}"
(cd "$root" && cargo build --release -j "${CARGO_BUILD_JOBS:-4}" --target wasm32-unknown-unknown -p rashevsky-wasm)
mkdir -p "$here/dist/fonts"
wasm-bindgen --target web --no-typescript --out-dir "$here/dist" "$target/wasm32-unknown-unknown/release/rashevsky_wasm.wasm"
if command -v wasm-opt > /dev/null; then
  wasm-opt -O3 --enable-bulk-memory --enable-nontrapping-float-to-int "$here/dist/rashevsky_wasm_bg.wasm" -o "$here/dist/rashevsky_wasm_bg.wasm"
fi
cp "$here/index.html" "$here/demo.js" "$here/demo.css" "$here/social.png" "$here/dist/"
# The MIT license of iCH360 asks for its notice in every copy, and the demo is one.
cp "$root/models/ich360/NOTICE" "$here/dist/NOTICE-iCH360.txt"
cp "$root/LICENSE" "$here/dist/LICENSE.txt"
cp "$root/THIRD-PARTY-NOTICES.md" "$here/dist/THIRD-PARTY-NOTICES.txt"
cp "$root/desktop/ui/js/fluxmap.js" "$root/desktop/ui/js/core.js" "$here/dist/"
cp "$root/desktop/ui/fonts/public-sans-wght-normal.woff2" "$root/desktop/ui/fonts/public-sans-OFL.txt" "$here/dist/fonts/"
ls -la "$here/dist"
