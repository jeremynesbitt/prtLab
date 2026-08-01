#!/usr/bin/env python3
# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Verify that every MATLAB source file has an explicit Python port status."""

from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PYTHON_ROOT.parent
MANIFEST = Path(__file__).with_name("matlab_port_manifest.json")
MATLAB_ROOT = REPOSITORY_ROOT / "prtLab"


def resolve_symbol(path: str):
    module_name, attribute = path.rsplit(".", 1)
    return getattr(importlib.import_module(module_name), attribute)


def audit() -> list[str]:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    categories = ("parity", "internal", "deferred", "legacy")
    classified: dict[str, str] = {}
    errors: list[str] = []
    for category in categories:
        names = (
            manifest[category]
            if isinstance(manifest[category], list)
            else manifest[category].keys()
        )
        for name in names:
            if name in classified:
                errors.append(
                    f"{name} is classified as both {classified[name]} and {category}"
                )
            classified[name] = category

    matlab_files = {path.name for path in MATLAB_ROOT.glob("*.m")}
    for name in sorted(matlab_files - classified.keys()):
        errors.append(f"{name} is not classified")
    for name in sorted(classified.keys() - matlab_files):
        errors.append(f"{name} is classified but does not exist")
    for matlab_name, python_symbol in manifest["parity"].items():
        try:
            resolve_symbol(python_symbol)
        except (AttributeError, ImportError) as error:
            errors.append(f"{matlab_name}: cannot import {python_symbol}: {error}")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--list-deferred",
        action="store_true",
        help="print MATLAB APIs intentionally deferred from the Python port",
    )
    arguments = parser.parse_args()
    errors = audit()
    if errors:
        raise SystemExit("API parity audit failed:\n  " + "\n  ".join(errors))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    print(
        "API parity manifest is complete: "
        f"{len(manifest['parity'])} parity, "
        f"{len(manifest['internal'])} internal, "
        f"{len(manifest['deferred'])} deferred, "
        f"{len(manifest['legacy'])} legacy."
    )
    if arguments.list_deferred and manifest["deferred"]:
        print("Deferred MATLAB APIs:")
        for name in manifest["deferred"]:
            print(f"  {name}")


if __name__ == "__main__":
    main()
