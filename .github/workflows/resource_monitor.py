#!/usr/bin/env python3
"""Print compact one-line disk/memory snapshots until the surrounding step ends.

Usage - add ONE line at the top of the workflow step whose resources you
want to watch (before any `cd`), leaving the rest of the step untouched:

    run: |
      python3 .github/workflows/resource_monitor.py &

      ... the actual build commands, unchanged ...

Output - a single colored line every --interval seconds (default 20),
interleaved with the build log but visually distinct and greppable:

    [RESOURCE MONITOR] 14:32:10 UTC | disk: 38G/145G (27%) - 107G free | mem: 1.3G/15G (8%) - 14G avail

The line color escalates with disk usage: cyan (normal) -> yellow (>= 75%)
-> red (>= 90%). Grep the raw log for "RESOURCE MONITOR" to extract the
full time series.

Annotations - the first time disk usage crosses a threshold, the script
emits a GitHub workflow-command annotation, which renders in the
Annotations panel on the run's summary page (separate from the logs):

    disk >= 75%  ->  ::warning::   (once per run)
    disk >= 90%  ->  ::error::     (once per run)
    mem available < 10%  ->  ::warning::  (once per run)

Annotations are processed live as the log streams, so they survive a runner
shutdown. They fire on threshold *crossings* only (not every sample),
because GitHub caps annotations at ~10 per type per step.
"""

import argparse
import os
import shutil
import sys
import time

CYAN = "\033[1;36m"
YELLOW = "\033[1;33m"
RED = "\033[1;31m"
RESET = "\033[0m"

DISK_WARN_PCT = 75
DISK_ERROR_PCT = 90
MEM_AVAIL_WARN_PCT = 10


def human(nbytes):
    """df-style short sizes: 38G, 1.3G, 512M, 107G."""
    for unit in ("B", "K", "M", "G", "T"):
        if nbytes < 1024 or unit == "T":
            if unit != "B" and nbytes < 10:
                return "{:.1f}{}".format(nbytes, unit)
            return "{:.0f}{}".format(nbytes, unit)
        nbytes /= 1024.0


def mem_info():
    """Return (total_bytes, available_bytes) from /proc/meminfo."""
    total = avail = 0
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemTotal:"):
                total = int(line.split()[1]) * 1024
            elif line.startswith("MemAvailable:"):
                avail = int(line.split()[1]) * 1024
    return total, avail


def emit(text):
    print(text, flush=True)


def snapshot(state):
    du = shutil.disk_usage("/")
    disk_pct = int(round(du.used * 100.0 / du.total))
    mem_total, mem_avail = mem_info()
    mem_used = mem_total - mem_avail
    mem_pct = int(round(mem_used * 100.0 / mem_total)) if mem_total else 0
    mem_avail_pct = 100 - mem_pct
    ts = time.strftime("%H:%M:%S", time.gmtime())

    color = RED if disk_pct >= DISK_ERROR_PCT else \
        YELLOW if disk_pct >= DISK_WARN_PCT else CYAN
    emit("{}[RESOURCE MONITOR] {} UTC \u2502 disk: {}/{} ({}%) \u2014 {} free"
         " \u2502 mem: {}/{} ({}%) \u2014 {} avail{}".format(
             color, ts,
             human(du.used), human(du.total), disk_pct, human(du.free),
             human(mem_used), human(mem_total), mem_pct, human(mem_avail),
             RESET))

    # Threshold-crossing annotations (plain lines, no ANSI codes: workflow
    # commands must start at column 0 with '::' to be parsed by the runner).
    if disk_pct >= DISK_ERROR_PCT and not state["disk_error"]:
        state["disk_error"] = True
        state["disk_warn"] = True
        emit("::error title=Runner disk almost full::disk / reached {}% "
             "({} free) at {} UTC \u2014 the runner may be killed if the disk "
             "fills completely (e.g. during docker layer export)".format(
                 disk_pct, human(du.free), ts))
    elif disk_pct >= DISK_WARN_PCT and not state["disk_warn"]:
        state["disk_warn"] = True
        emit("::warning title=Runner disk usage high::disk / reached {}% "
             "({} free) at {} UTC".format(disk_pct, human(du.free), ts))

    if mem_avail_pct < MEM_AVAIL_WARN_PCT and not state["mem_warn"]:
        state["mem_warn"] = True
        emit("::warning title=Runner memory low::only {} of memory available "
             "({}% used) at {} UTC".format(human(mem_avail), mem_pct, ts))


def main():
    parser = argparse.ArgumentParser(
        description="Background resource monitor for GitHub Actions steps.")
    parser.add_argument("--interval", type=float, default=20,
                        help="seconds between snapshots (default: 20)")
    args = parser.parse_args()

    parent_pid = os.getppid()
    state = {"disk_warn": False, "disk_error": False, "mem_warn": False}

    try:
        emit("{}[RESOURCE MONITOR] started \u2014 sampling every {:g} s{}".format(
            CYAN, args.interval, RESET))
        snapshot(state)
        while True:
            time.sleep(args.interval)
            if os.getppid() != parent_pid:
                # The step's shell exited and we were re-parented: stop quietly.
                return
            snapshot(state)
    except (BrokenPipeError, KeyboardInterrupt):
        # Runner closed the log pipe (step over) or the job was cancelled.
        # Suppress the BrokenPipeError noise Python emits on interpreter
        # shutdown when stdout is already gone.
        try:
            sys.stdout.close()
        except Exception:
            pass


if __name__ == "__main__":
    main()
