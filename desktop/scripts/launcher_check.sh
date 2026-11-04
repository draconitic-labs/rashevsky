#!/usr/bin/env bash
# Starts the installed entry with `gtk-launch` inside a virtual KWin session, a real window manager on a
# private display, away from the desktop. Reports the window's class, name and icon, closes it through
# the window manager, and leaves nothing running. There is no session bus of its own, because the RustDesk
# service restarts itself when a new dbus-daemon appears. A virtual KWin still makes that service restart
# once, which drops any remote session, so do not run this while connected through RustDesk.
set -euo pipefail
out="${1:-$HOME/.cache/rashevsky-work/launcher}"
mkdir -p "$out"
if [[ "${RASHEVSKY_LAUNCHER_INNER:-}" != "1" ]]; then
  for entry in "$HOME/.local/share/applications/rashevsky.desktop" "$(xdg-user-dir DESKTOP)/Rashevsky.desktop"; do
    desktop-file-validate "$entry" && echo "valid: $entry, mode $(stat -c %a "$entry")"
  done
  gio info -a metadata::trusted "$(xdg-user-dir DESKTOP)/Rashevsky.desktop" | tail -1
  export RASHEVSKY_LAUNCHER_INNER=1
  env -u WAYLAND_DISPLAY -u DISPLAY GDK_BACKEND=x11 DBUS_SESSION_BUS_ADDRESS=unix:path=/nonexistent \
    timeout 60 kwin_wayland --virtual --xwayland --no-lockscreen --width 1600 --height 1000 \
    --exit-with-session "$0 $out" > "$out/kwin.log" 2>&1 || true
  grep -v -E 'a11y|AT-SPI|^$' "$out/report.txt"
  pgrep -x rashevsky-deskt > /dev/null && echo "process left running" || echo "no process left"
  exit 0
fi
report="$out/report.txt"
: > "$report"
for _ in $(seq 1 60); do xdpyinfo > /dev/null 2>&1 && break; sleep 0.5; done
gtk-launch rashevsky > "$out/launch.log" 2>&1
echo "gtk-launch exit $?" >> "$report"
window=""
for _ in $(seq 1 80); do
  window="$(xdotool search --onlyvisible --name '^Rashevsky$' 2> /dev/null | head -1)" && [[ -n "$window" ]] && break
  sleep 0.25
done
echo "window $window" >> "$report"
sleep 5
xprop -id "$window" WM_CLASS WM_NAME >> "$report" 2>&1
icon_values="$(xprop -id "$window" -notype 32c _NET_WM_ICON | tr -cd ',' | wc -c)"
echo "_NET_WM_ICON values $icon_values" >> "$report"
wmctrl -lx >> "$report" 2>&1
import -window "$window" "$out/window.png" 2>> "$report" || true
wmctrl -c Rashevsky >> "$report" 2>&1
for _ in $(seq 1 40); do pgrep -x rashevsky-deskt > /dev/null || break; sleep 0.25; done
if pgrep -x rashevsky-deskt > /dev/null; then echo "after close: running" >> "$report"; else echo "after close: exited" >> "$report"; fi
