#!/usr/bin/env python3
# BRL-CAD
#
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# version 2.1 as published by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this file; see the file named COPYING for more
# information.
#
"""Run configured MGED tests with a substitute frontend in an isolated workspace."""

import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import tempfile
import time


def mged_test(test):
    for argument in test.get("command", []):
        if re.search(r"(?:^|/|=)mged(?:\.exe)?$", argument):
            return True
        script = Path(argument)
        if script.suffix == ".sh" and script.is_file():
            if re.search(r"ensearch[^\n]*mged|\$\{?MGED\b", script.read_text()):
                return True
    return False


def test_properties(test):
    return {entry["name"]: entry["value"] for entry in test["properties"]}


def fixture_tests(inventory, selected):
    setups = {}
    for test in inventory:
        for fixture in test_properties(test).get("FIXTURES_SETUP", []):
            setups.setdefault(fixture, []).append(test)
    tests = {test["name"]: test for test in selected}
    required = set()
    while True:
        for test in list(tests.values()):
            required.update(test_properties(test).get("FIXTURES_REQUIRED", []))
        additions = {test["name"]: test for fixture in required
                     for test in setups.get(fixture, []) if test["name"] not in tests}
        if not additions:
            break
        tests.update(additions)
    for test in tests.values():
        required.update(test_properties(test).get("FIXTURES_SETUP", []))
    cleanup = [test for test in inventory
               if required.intersection(test_properties(test).get("FIXTURES_CLEANUP", []))]
    return list(tests.values()), cleanup


def run_tests(tests, build, workspace, wrapper, jobs, timeout):
    pending = {test["name"]: test for test in tests}
    properties_by_name = {test["name"]: test_properties(test) for test in tests}
    dependencies = {}
    required_setups = {}
    for test in tests:
        properties = properties_by_name[test["name"]]
        required = set(properties.get("FIXTURES_REQUIRED", []))
        cleanup = set(properties.get("FIXTURES_CLEANUP", []))
        dependencies[test["name"]] = set(properties.get("DEPENDS", [])).intersection(pending)
        required_setups[test["name"]] = set()
        for other_name, other_properties in properties_by_name.items():
            if other_name == test["name"]:
                continue
            if required.intersection(other_properties.get("FIXTURES_SETUP", [])):
                dependencies[test["name"]].add(other_name)
                required_setups[test["name"]].add(other_name)
            if cleanup.intersection(
                    other_properties.get("FIXTURES_SETUP", []) + other_properties.get("FIXTURES_REQUIRED", [])):
                dependencies[test["name"]].add(other_name)
    completed = {}
    running = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        while pending or running:
            made_progress = False
            for name, test in list(pending.items()):
                if len(running) >= jobs:
                    break
                if dependencies[name].issubset(completed):
                    failed_setups = [setup for setup in required_setups[name] if completed[setup]["status"]]
                    if failed_setups:
                        reason = "Fixture setup failed: " + ", ".join(sorted(failed_setups))
                        print(f"FAIL   {name} (not run: {reason})", flush=True)
                        (workspace / (name + ".txt")).write_text(reason + "\n")
                        completed[name] = {"name": name, "outcome": "fail", "status": "not run",
                                           "invocations": 0, "seconds": 0, "reason": reason}
                        del pending[name]
                        made_progress = True
                        continue
                    future = pool.submit(run_test, test, build, workspace, wrapper, timeout)
                    running[future] = name
                    del pending[name]
            if not running:
                if not pending:
                    break
                if made_progress:
                    continue
                raise RuntimeError("Regression fixture dependencies contain a cycle")
            done, _ = concurrent.futures.wait(running, return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                name = running.pop(future)
                completed[name] = future.result()
    return [completed[test["name"]] for test in tests]


def run_test(test, build, workspace, wrapper, timeout):
    name = test["name"]
    properties = test_properties(test)
    original_working = Path(properties.get("WORKING_DIRECTORY", build))
    try:
        relative = original_working.relative_to(build)
    except ValueError:
        relative = Path(name)
    working = workspace / relative
    working.mkdir(parents=True, exist_ok=True)
    trace = workspace / (name + ".trace")
    original_mged = str(build / "bin" / wrapper.name)
    command = [argument.replace(original_mged, str(wrapper)) for argument in test["command"]]
    if properties.get("FIXTURES_CLEANUP"):
        command = [argument.replace(str(build) + os.sep, str(workspace) + os.sep)
                   for argument in command]
    for index, argument in enumerate(command):
        if argument.startswith("-DTEST_DIR="):
            test_directory = Path(argument.split("=", 1)[1])
            try:
                relative_directory = test_directory.relative_to(build)
            except ValueError:
                continue
            command[index] = "-DTEST_DIR=" + str(workspace / relative_directory)

    # These configured scripts embed output directories. Keep their input
    # scripts in the build tree and move only generated files into the audit.
    if name.startswith("regress-mged") and Path(command[-1]).suffix == ".cmake" and not properties.get("FIXTURES_CLEANUP"):
        configured = Path(command[-1])
        contents = configured.read_text()
        contents = contents.replace('set(CBDIR "' + str(original_working),
                                    'set(CBDIR "' + str(working))
        contents = contents.replace('"${CBDIR}/regress_${TNAME}.mged"',
                                    '"' + str(original_working) + '/regress_${TNAME}.mged"')
        contents = re.sub(r'set\(LOGFILE "[^"\n]*"\)',
                          lambda _: 'set(LOGFILE "' + str(working / (name + ".log")) + '")', contents)
        configured_copy = workspace / (name + ".cmake")
        configured_copy.write_text(contents)
        command[-1] = str(configured_copy)

    environment = os.environ.copy()
    environment["FLTKGED_AUDIT_TRACE"] = str(trace)
    environment["PATH"] = str(wrapper.parent) + os.pathsep + environment.get("PATH", "")
    for entry in properties.get("ENVIRONMENT", []):
        key, value = entry.split("=", 1)
        environment[key] = value
    for key, directory in (("BU_DIR_CACHE", "cache"), ("LIBRT_CACHE", "rtcache")):
        cache = working / (name + "-" + directory)
        cache.mkdir(exist_ok=True)
        environment[key] = str(cache)

    start = time.monotonic()
    process = subprocess.Popen(command, cwd=working, env=environment, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, start_new_session=os.name == "posix")
    try:
        output, _ = process.communicate(timeout=min(properties.get("TIMEOUT", timeout), timeout))
        status = process.returncode
    except subprocess.TimeoutExpired:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGKILL)
        else:
            process.kill()
        output, _ = process.communicate()
        status = 124
    invocations = trace.read_text().count("invoked") if trace.exists() else 0
    log = working / (name + ".log")
    if log.is_file():
        output += "\n" + log.read_text(errors="replace")
    (workspace / (name + ".txt")).write_text(output)
    outcome = "unused" if not invocations and status == 0 else "pass" if status == 0 else "fail"
    label = "FIXTURE" if properties.get("FIXTURES_CLEANUP") else outcome.upper()
    print(f"{label:7} {name} (status {status}, {invocations} invocations)", flush=True)
    return {"name": name, "outcome": outcome, "status": status,
            "invocations": invocations, "seconds": round(time.monotonic() - start, 2)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_directory", type=Path)
    parser.add_argument("--frontend", default="fltkged", help="Build executable name or absolute path")
    parser.add_argument("--output-directory", type=Path)
    parser.add_argument("--include", default=".", help="Regular expression selecting test names")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=300)
    options = parser.parse_args()
    build = options.build_directory.resolve()
    frontend = Path(options.frontend)
    if not frontend.is_absolute():
        frontend = build / "bin" / frontend
    if not frontend.is_file():
        parser.error(f"Frontend does not exist: {frontend}")
    if options.jobs < 1 or options.timeout <= 0:
        parser.error("Jobs and timeout must be positive")
    selection = re.compile(options.include)
    inventory = subprocess.run(["ctest", "--test-dir", str(build), "--show-only=json-v1"],
                               check=True, capture_output=True, text=True)
    inventory = json.loads(inventory.stdout)["tests"]
    tests = [test for test in inventory
             if mged_test(test) and selection.search(test["name"])]
    tests, cleanup = fixture_tests(inventory, tests)
    for test in tests:
        # ensearch checks the source tree's bin directory before the working
        # tree. Refuse a setup where that could bypass the substitute frontend.
        for argument in test["command"]:
            if argument.startswith("-"):
                continue
            source = Path(argument)
            if source.is_dir() and (source / "bin" / "mged").is_file() and source != build:
                parser.error("An in-source MGED executable would bypass the audit wrapper")
    parent = options.output_directory
    if parent:
        parent.mkdir(parents=True, exist_ok=True)
    workspace = Path(tempfile.mkdtemp(prefix="fltkged-regress-", dir=parent))
    binaries = workspace / "bin"
    binaries.mkdir()
    executable_suffix = frontend.suffix if frontend.suffix == ".exe" else ""
    wrapper = binaries / ("mged" + executable_suffix)
    if os.name != "posix":
        parser.error("The shell-based regression audit currently requires a POSIX host")
    for original in (build / "bin").iterdir():
        if original.name != wrapper.name:
            (binaries / original.name).symlink_to(original)
    for directory in ("lib", "share", "include"):
        if (build / directory).exists():
            (workspace / directory).symlink_to(build / directory, target_is_directory=True)
    wrapper.write_text('#!/bin/sh\necho invoked >> "$FLTKGED_AUDIT_TRACE"\nexec ' +
                       shlex.quote(str(frontend)) + ' "$@"\n')
    wrapper.chmod(0o755)
    print(f"Audit workspace: {workspace}", flush=True)
    outcomes = run_tests(tests + cleanup, build, workspace, wrapper, options.jobs, options.timeout)
    results = outcomes[:len(tests)]
    counts = {outcome: sum(result["outcome"] == outcome for result in results)
              for outcome in ("pass", "fail", "unused")}
    report = {"frontend": frontend.name, "counts": counts, "tests": results,
              "fixture_cleanup": outcomes[len(tests):]}
    (workspace / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Passed {counts['pass']}; failed {counts['fail']}; did not invoke MGED {counts['unused']}")
    print(f"Results: {workspace / 'results.json'}")
    return int(counts["fail"] != 0 or any(result["status"] for result in report["fixture_cleanup"]))


if __name__ == "__main__":
    raise SystemExit(main())
