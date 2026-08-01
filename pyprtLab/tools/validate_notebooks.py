#!/usr/bin/env python3
# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Validate or execute every Python documentation notebook."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_DIR = ROOT / "docs" / "source"


def code_cells(path: Path):
    notebook = json.loads(path.read_text(encoding="utf-8"))
    if notebook.get("nbformat") != 4:
        raise ValueError(f"{path}: expected notebook format 4")
    for index, cell in enumerate(notebook.get("cells", []), start=1):
        if cell.get("cell_type") == "code":
            yield index, "".join(cell.get("source", []))


def validate(path: Path, execute: bool) -> None:
    scope = {"__name__": "__notebook__"}
    for index, source in code_cells(path):
        compiled = compile(source, f"{path.name}:cell-{index}", "exec")
        if execute:
            # This tool intentionally executes trusted, repository-owned notebooks.
            exec(compiled, scope)  # noqa: S102


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--execute", action="store_true",
        help="execute cells in addition to checking syntax",
    )
    arguments = parser.parse_args()
    if arguments.execute:
        import matplotlib

        matplotlib.use("Agg")
        warnings.filterwarnings(
            "ignore", message="FigureCanvasAgg is non-interactive"
        )
    notebooks = sorted(NOTEBOOK_DIR.glob("*.ipynb"))
    if not notebooks:
        raise SystemExit(f"No notebooks found in {NOTEBOOK_DIR}")
    for notebook in notebooks:
        validate(notebook, arguments.execute)
        print(f"{'Executed' if arguments.execute else 'Validated'} {notebook.name}")


if __name__ == "__main__":
    main()
