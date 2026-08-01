# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Load MATLAB reference fixtures with provenance checks."""

from __future__ import annotations

import json
from pathlib import Path

BASELINES = Path(__file__).parent / "baselines"
REQUIRED_METADATA = {
    "schemaVersion",
    "referenceImplementation",
    "referenceCommit",
    "generator",
    "absoluteTolerance",
}


def load_json_baseline(name: str) -> dict:
    data = json.loads((BASELINES / name).read_text(encoding="utf-8"))
    missing = REQUIRED_METADATA - data.keys()
    if missing:
        raise ValueError(f"{name} is missing baseline metadata: {sorted(missing)}")
    if len(data["referenceCommit"]) != 40:
        raise ValueError(f"{name} has an invalid reference commit")
    return data
