#!/usr/bin/env bash
# Renders the application icons from their SVG sources with Inkscape. Sizes of 32 and below use the
# simplified drawing, in which the beads are left out so the letter stays readable.
set -euo pipefail
icons="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../src-tauri/icons" && pwd)"
for size in 16 24 32; do
  inkscape "$icons/rashevsky-small.svg" --export-type=png -w "$size" -h "$size" -o "$icons/icon-$size.png" 2>/dev/null
done
for size in 48 64 128 256 512; do
  inkscape "$icons/rashevsky.svg" --export-type=png -w "$size" -h "$size" -o "$icons/icon-$size.png" 2>/dev/null
done
cp "$icons/icon-512.png" "$icons/icon.png"
