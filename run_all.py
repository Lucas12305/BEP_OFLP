#!/usr/bin/env python
# coding: utf-8

# # Running all of the instances in all of the settings
# This document runs every instance with every combination of method (LP or GS), precutting (off or on) and
# separator (off or on). Every run is a separate call of run_OFLP.py, so that a crash or a time out of one run
# does not stop the others. The output of each run is saved in a text file in the results folder, and after
# every run a line is added to a summary csv file, so the progress is never lost.
#
# Usage: python run_all.py [--time-limit SECONDS] [--hard-timeout SECONDS] [--dry-run] [instance ...]
#
# Examples:
#   python run_all.py --time-limit 600
#   python run_all.py --time-limit 600 ofl_2_100i ofl_3_100i
#   python run_all.py --dry-run

import argparse
import csv
import itertools
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# The folder of this file, so that the script can be called from anywhere
HERE = Path(__file__).resolve().parent
INSTANCE_DIR = HERE / "initialized_instances"
OUT_DIR = HERE / "results"

MODES = [0, 1]
# 0 is LP, 1 is GS
PRECUTTINGS = [0, 1]
# 0 is off, 1 is on
SEPERATORS = [0, 1]
# 0 is off, 1 is on

STATUS_RE = re.compile(r"^Status: (.+?), total time: ([\d.]+) s", re.MULTILINE)
OBJECTIVE_RE = re.compile(r"^Objective L \(squared distance\): (.+)$", re.MULTILINE)


def natural_key(name : str) -> list:
    '''
    This definition gives a key to sort names naturally, so that ofl_2_100i comes before ofl_10_100i.
    :param name: The name that is sorted.
    '''
    parts = re.split(r"(\d+)", name)
    key = []
    for p in parts:
        if p.isdigit():
            key.append(int(p))
        else:
            key.append(p)
    return key


def find_all_instances() -> list:
    '''
    This definition finds the names (without .gms) of all the instances in the initialized_instances folder.
    '''
    names = []
    for path in INSTANCE_DIR.glob("*.gms"):
        names.append(path.stem)
    names.sort(key=natural_key)
    return names


def read_summary(text : str) -> tuple:
    '''
    This definition reads the status, the total time and the objective from the output of run_OFLP.py.
    When it is not in the output (for example when the run crashed) these are empty.
    :param text: The output of the run.
    '''
    status, solve_time, objective = "", "", ""
    found = STATUS_RE.search(text)
    if found:
        status, solve_time = found.group(1), found.group(2)
    found = OBJECTIVE_RE.search(text)
    if found:
        objective = found.group(1)
    return status, solve_time, objective


def run_one(name : str, mode : int, precutting : int, seperator : int, time_limit : float, hard_timeout : float, plot : bool, stamp : str) -> dict:
    '''
    This definition runs one instance with one setting, the output is saved in a text file in the results folder.
    :param name: The name of the instance, for example ofl_2_100i.
    :param mode: 0 is LP, 1 is GS.
    :param precutting: 0 is off, 1 is on.
    :param seperator: 0 is off, 1 is on.
    :param time_limit: The time limit for SCIP in seconds, when it is None the default of run_OFLP.py is used.
    :param hard_timeout: The time after which the run is killed in seconds, when it is None there is no limit.
    :param plot: If this is True the plot of the solution is saved next to the text file.
    :param stamp: The time stamp that is used in the names of the files.
    :return: A dictionary with the result of this run.
    '''
    base = OUT_DIR / f"{name}_mode{mode}_precut{precutting}_sep{seperator}_{stamp}"
    command = [sys.executable, "-u", str(HERE / "run_OFLP.py"), name,
               "--mode", str(mode), "--precutting", str(precutting), "--seperator", str(seperator)]
    if time_limit is not None:
        command += ["--time-limit", str(time_limit)]
    if plot:
        command += ["--plot-file", f"{base}.png"]

    start = time.perf_counter()
    try:
        done = subprocess.run(command, cwd=HERE, capture_output=True, text=True, timeout=hard_timeout)
        output = done.stdout + done.stderr
        if done.returncode == 0:
            outcome = "ok"
        else:
            outcome = f"failed (exit code {done.returncode})"
    except subprocess.TimeoutExpired as e:
        # With a time out, the output until now is still saved
        output = ""
        for part in (e.stdout, e.stderr):
            if part is not None:
                output += part.decode() if isinstance(part, bytes) else part
        outcome = "killed (hard timeout)"
    wall_time = time.perf_counter() - start

    Path(f"{base}.txt").write_text(output)
    status, solve_time, objective = read_summary(output)
    return {"instance": name, "mode": mode, "precutting": precutting, "seperator": seperator,
            "outcome": outcome, "scip_status": status, "scip_time_s": solve_time,
            "wall_time_s": f"{wall_time:.1f}", "objective_L": objective, "output_file": f"{base.name}.txt"}


def main() -> None:
    '''
    This reads the command line and runs all of the instances in all of the settings.
    '''
    parser = argparse.ArgumentParser(description="Run all instances in all settings (LP/GS, precutting, seperator).")
    parser.add_argument("instances", nargs="*", help="The names of the instances, when none are given all instances are run")
    parser.add_argument("--time-limit", type=float, default=600, help="The time limit for SCIP in seconds per run (default 600)")
    parser.add_argument("--hard-timeout", type=float, default=None, help="Kill a run after this many seconds, this also covers the steps before SCIP starts (default no limit)")
    parser.add_argument("--plot", action="store_true", help="Also save a plot of the solution of every run")
    parser.add_argument("--dry-run", action="store_true", help="Only print the runs that would be done")
    args = parser.parse_args()

    names = args.instances
    if len(names) == 0:
        names = find_all_instances()
    for i in range(len(names)):
        names[i] = Path(names[i]).stem
        if not (INSTANCE_DIR / f"{names[i]}.gms").exists():
            print(f"Instance {names[i]} was not found in {INSTANCE_DIR}.", file=sys.stderr)
            sys.exit(1)

    settings = list(itertools.product(MODES, PRECUTTINGS, SEPERATORS))
    total = len(names) * len(settings)
    print(f"{len(names)} instances x {len(settings)} settings = {total} runs, time limit {args.time_limit} s per run", flush=True)

    if args.dry_run:
        for name in names:
            for mode, precutting, seperator in settings:
                print(f"{name}: mode {mode}, precutting {precutting}, seperator {seperator}")
        return None

    OUT_DIR.mkdir(exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    summary_path = OUT_DIR / f"summary_{stamp}.csv"
    fields = ["instance", "mode", "precutting", "seperator", "outcome", "scip_status",
              "scip_time_s", "wall_time_s", "objective_L", "output_file"]

    count = 0
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for name in names:
            for mode, precutting, seperator in settings:
                count += 1
                print(f"[{count}/{total}] {name}: mode {mode}, precutting {precutting}, seperator {seperator} ...", end=" ", flush=True)
                result = run_one(name, mode, precutting, seperator, args.time_limit, args.hard_timeout, args.plot, stamp)
                print(f"{result['outcome']} in {result['wall_time_s']} s", flush=True)
                writer.writerow(result)
                f.flush()
    print(f"Done, the summary is saved in {summary_path}")
    return None


if __name__ == "__main__":
    main()