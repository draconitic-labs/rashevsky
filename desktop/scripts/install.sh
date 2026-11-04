#!/usr/bin/env bash
# Builds the release binary and installs it for this user, with its icons, its application entry and
# a shortcut on the desktop. Writes only Rashevsky's own files. Nothing here needs root.
set -euo pipefail

source_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
tauri_dir="$source_root/src-tauri"
target_dir="${CARGO_TARGET_DIR:-$tauri_dir/target}"

if [[ "${1:-}" != "--no-build" ]]; then
  (cd "$tauri_dir" && CARGO_TARGET_DIR="$target_dir" cargo build --release)
fi

bin_dir="$HOME/.local/bin"
apps_dir="$HOME/.local/share/applications"
icon_root="$HOME/.local/share/icons/hicolor"
desktop_dir="$(xdg-user-dir DESKTOP 2>/dev/null || echo "$HOME/Desktop")"
binary="$bin_dir/rashevsky-desktop"

mkdir -p "$bin_dir" "$apps_dir" "$desktop_dir"
# Replacing rather than overwriting keeps a running copy of the application intact.
install -m 755 "$target_dir/release/rashevsky-desktop" "$binary.new"
mv -f "$binary.new" "$binary"

for size in 16 24 32 48 64 128 256 512; do
  mkdir -p "$icon_root/${size}x${size}/apps"
  install -m 644 "$tauri_dir/icons/icon-$size.png" "$icon_root/${size}x${size}/apps/rashevsky.png"
done
mkdir -p "$icon_root/scalable/apps"
install -m 644 "$tauri_dir/icons/rashevsky.svg" "$icon_root/scalable/apps/rashevsky.svg"

entry="$(mktemp --tmpdir="${TMPDIR:-/tmp}" --suffix=.desktop rashevsky-entry.XXXXXX)"
trap 'unlink "$entry" 2>/dev/null || true' EXIT
cat > "$entry" <<DESKTOP
[Desktop Entry]
Type=Application
Version=1.5
Name=Rashevsky
GenericName=Metabolic network modelling
Comment=Watch flux run through metabolism, cut a reaction, and see the cell reroute. Drop an SBML file here to open it.
Exec=$binary %f
TryExec=$binary
Icon=rashevsky
Terminal=false
Categories=Science;Biology;Education;
Keywords=metabolism;flux balance;FBA;SBML;knockout;kinetics;biochemistry;systems biology;COBRA;
MimeType=application/sbml+xml;application/xml;text/xml;application/gzip;
InitialPreference=1
StartupNotify=true
StartupWMClass=rashevsky-desktop
Actions=screen;kinetics;atlas;checks;

[Desktop Action screen]
Name=Run a knockout screen
Exec=$binary --view=screen
Icon=rashevsky

[Desktop Action kinetics]
Name=Simulate yeast glycolysis
Exec=$binary --view=kinetics
Icon=rashevsky

[Desktop Action atlas]
Name=Fluxes on the human atlas
Exec=$binary --view=atlas
Icon=rashevsky

[Desktop Action checks]
Name=Show validation results
Exec=$binary --view=checks
Icon=rashevsky
DESKTOP
desktop-file-validate "$entry"

install -m 644 "$entry" "$apps_dir/rashevsky.desktop"
install -m 755 "$entry" "$desktop_dir/Rashevsky.desktop"
# GNOME and KDE both respect this flag; without it they ask before running a desktop file.
gio set "$desktop_dir/Rashevsky.desktop" metadata::trusted true 2>/dev/null || true
update-desktop-database "$apps_dir" 2>/dev/null || true
gtk-update-icon-cache -q -t "$icon_root" 2>/dev/null || true
echo "Installed $binary and $desktop_dir/Rashevsky.desktop"
