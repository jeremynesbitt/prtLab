#!/usr/bin/env python3
# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Report MATLAB changes since the last acknowledged Python sync commit."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PYTHON_ROOT.parent
MARKER = PYTHON_ROOT / "MATLAB_SYNC_COMMIT"
MATLAB_PATHS = ("prtLab", "examples", "tests", "docs/source")
IGNORED_WORKTREE_SUFFIXES = (".asv", ".DS_Store")
MANIFEST = Path(__file__).with_name("matlab_port_manifest.json")


def git(*arguments: str) -> str:
    result = subprocess.run(
        ("git", *arguments),
        cwd=REPOSITORY_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def relevant_status(status: str) -> str:
    return "\n".join(
        line for line in status.splitlines()
        if not line.rstrip().endswith(IGNORED_WORKTREE_SUFFIXES)
    )


def port_status(path: str, manifest: dict) -> str:
    name = Path(path).name
    for category in ("parity", "internal", "deferred", "legacy"):
        if name in manifest[category]:
            return category
    return "unclassified"


def annotate_changes(changes: str, manifest: dict) -> str:
    annotated = []
    for line in changes.splitlines():
        path = line.split()[-1]
        status = port_status(path, manifest) if path.startswith("prtLab/") else "support"
        annotated.append(f"{line}  [{status}]")
    return "\n".join(annotated)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mark-current",
        action="store_true",
        help="record the current HEAD after Python parity has been reviewed",
    )
    arguments = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    head = git("rev-parse", "HEAD")
    if arguments.mark_current:
        MARKER.write_text(head + "\n", encoding="utf-8")
        print(f"Recorded MATLAB sync commit {head}")
        return

    baseline = MARKER.read_text(encoding="utf-8").strip()
    if baseline == "UNINITIALIZED":
        raise SystemExit(
            "MATLAB_SYNC_COMMIT is uninitialized. Run this tool with "
            "--mark-current after reviewing the initial Python port."
        )
    committed = git(
        "diff", "--name-status", f"{baseline}..HEAD", "--", *MATLAB_PATHS
    )
    working = relevant_status(git("status", "--short", "--", *MATLAB_PATHS))
    if not committed and not working:
        print(f"Python is synchronized with tracked MATLAB changes through {head}.")
        return
    if committed:
        print("Committed MATLAB changes since the Python sync point:")
        print(annotate_changes(committed, manifest))
    if working:
        print("Uncommitted MATLAB changes:")
        print(annotate_changes(working, manifest))
    raise SystemExit(1)


if __name__ == "__main__":
    main()
