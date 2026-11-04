#!/usr/bin/env python3
"""Drives the real Tauri window on a private X display with real keyboard and pointer input.

Everything runs on its own Xvfb display with no session bus, so the desktop is not touched. The application
writes frame rates to its log when RASHEVSKY_QA is set.
"""
from __future__ import annotations

import os
import pathlib
import select
import signal
import subprocess
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
BINARY = pathlib.Path(os.environ.get("CARGO_TARGET_DIR", ROOT / "src-tauri" / "target")) / "release" / "rashevsky-desktop"
SCRATCH = pathlib.Path(os.environ.get("XDG_CACHE_HOME", pathlib.Path.home() / ".cache")) / "rashevsky-qa"


class Session:
    def __init__(self, size=(1560, 980), screen=(1760, 1200), args=(), binary=BINARY, title="^Rashevsky$"):
        SCRATCH.mkdir(parents=True, exist_ok=True)
        (SCRATCH / "tmp").mkdir(exist_ok=True)
        self.log_path = SCRATCH / f"app-{time.time_ns()}.log"
        self.log = self.log_path.open("w")
        read_fd, write_fd = os.pipe()
        self.xvfb = subprocess.Popen(
            ["Xvfb", "-displayfd", str(write_fd), "-screen", "0", f"{screen[0]}x{screen[1]}x24", "-nolisten", "tcp"],
            pass_fds=(write_fd,), stdout=self.log, stderr=self.log, start_new_session=True)
        os.close(write_fd)
        # Xvfb writes the number and the newline separately and dies if the pipe has closed in between.
        number = b""
        while not number.endswith(b"\n"):
            if not select.select([read_fd], [], [], 15)[0]:
                raise RuntimeError("the private display did not start")
            chunk = os.read(read_fd, 100)
            if not chunk:
                raise RuntimeError(f"the private display exited; see {self.log_path}")
            number += chunk
        self.display = ":" + number.decode().strip()
        os.close(read_fd)
        self.env = os.environ.copy()
        self.env.update(DISPLAY=self.display, GDK_BACKEND="x11", WEBKIT_DISABLE_COMPOSITING_MODE="1",
                        WEBKIT_DISABLE_DMABUF_RENDERER="1", TMPDIR=str(SCRATCH / "tmp"), RASHEVSKY_QA="1")
        self.env.pop("WAYLAND_DISPLAY", None)
        # No session bus at all. A private `dbus-run-session` is not an option here: the RustDesk service
        # follows the newest dbus-daemon and restarts its server and tray whenever one appears or goes.
        self.env.update(DBUS_SESSION_BUS_ADDRESS="unix:path=/nonexistent", NO_AT_BRIDGE="1", GTK_A11Y="none")
        self.app = subprocess.Popen([str(binary), *args], env=self.env, stdout=self.log,
                                    stderr=self.log, start_new_session=True)
        self.title = title
        self.window = self.wait_window()
        self.xdo("windowsize", self.window, str(size[0]), str(size[1]))
        self.xdo("windowmove", self.window, "0", "0")
        time.sleep(1.2)

    def wait_window(self):
        for _ in range(240):
            found = subprocess.run(["xdotool", "search", "--onlyvisible", "--name", self.title], env=self.env, capture_output=True, text=True)
            if found.stdout.strip():
                return found.stdout.split()[0]
            if self.app.poll() is not None:
                raise RuntimeError(f"the application exited {self.app.returncode}; see {self.log_path}")
            time.sleep(0.25)
        raise RuntimeError("no window appeared")

    def xdo(self, *args):
        subprocess.run(["xdotool", *args], env=self.env, check=True)

    def key(self, *keys, settle=0.5):
        self.xdo("windowfocus", self.window)
        for k in keys:
            self.xdo("key", "--window", self.window, k)
        time.sleep(settle)

    def type(self, text, settle=0.5):
        self.xdo("type", "--window", self.window, text)
        time.sleep(settle)

    def move(self, x, y, settle=0.3):
        self.xdo("mousemove", "--window", self.window, str(x), str(y))
        time.sleep(settle)

    def click(self, x, y, settle=0.6, button=1, repeat=1):
        self.move(x, y, 0.15)
        self.xdo("click", "--repeat", str(repeat), "--delay", "90", str(button))
        time.sleep(settle)

    def shot(self, path):
        path = pathlib.Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["import", "-window", self.window, str(path)], env=self.env, check=True)
        return path

    def read_log(self):
        self.log.flush()
        return self.log_path.read_text(errors="replace")

    def close(self):
        for proc in (self.app, self.xvfb):
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        for proc in (self.app, self.xvfb):
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
        self.log.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
