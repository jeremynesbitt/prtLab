# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import json
from pathlib import Path

import pytest

NOTEBOOKS = sorted((Path(__file__).parents[1] / "docs" / "source").glob("*.ipynb"))
REPOSITORY_ROOT = Path(__file__).parents[2]
PAIRED_NOTEBOOKS = {
    "QuarterWavePlateAnalysis.ipynb": "QuarterWavePlateAnalysis.ipynb",
    "CellPhoneLensSystemExample.ipynb": "CellPhoneLensSystemExample.ipynb",
    "CassegrainTelescopeExample.ipynb": "CassegrainTelescopeExample.ipynb",
    "FresnelRhombAnalysis.ipynb": "FresnelRhombExample.ipynb",
}
ALLOWED_MARKDOWN_REPLACEMENTS = {
    "This is implemented using the table type in matlab.": (
        "This is implemented using the OpticalSystem type in Python."
    ),
    "the trace*.m files": "the trace interface files",
}


@pytest.mark.parametrize("path", NOTEBOOKS, ids=lambda path: path.stem)
def test_notebook_code_cells_compile(path):
    notebook = json.loads(path.read_text(encoding="utf-8"))
    assert notebook["nbformat"] == 4
    assert notebook["metadata"]["kernelspec"]["language"] == "python"
    assert notebook["metadata"]["kernelspec"]["name"] == "prtlab-python"
    for index, cell in enumerate(notebook["cells"], start=1):
        if cell["cell_type"] == "code":
            compile("".join(cell["source"]), f"{path.name}:cell-{index}", "exec")


@pytest.mark.parametrize("python_name,matlab_name", PAIRED_NOTEBOOKS.items())
def test_translated_notebook_preserves_matlab_markdown(
    python_name, matlab_name
):
    python_notebook = json.loads(
        (NOTEBOOKS[0].parent / python_name).read_text(encoding="utf-8")
    )
    matlab_notebook = json.loads(
        (REPOSITORY_ROOT / "docs" / "source" / matlab_name).read_text(
            encoding="utf-8"
        )
    )
    assert [
        cell["cell_type"] for cell in python_notebook["cells"]
    ] == [cell["cell_type"] for cell in matlab_notebook["cells"]]
    python_markdown = [
        "".join(cell["source"])
        for cell in python_notebook["cells"]
        if cell["cell_type"] == "markdown"
    ]
    matlab_markdown = [
        "".join(cell["source"])
        for cell in matlab_notebook["cells"]
        if cell["cell_type"] == "markdown"
    ]
    for old, new in ALLOWED_MARKDOWN_REPLACEMENTS.items():
        matlab_markdown = [text.replace(old, new) for text in matlab_markdown]
    assert python_markdown == matlab_markdown
