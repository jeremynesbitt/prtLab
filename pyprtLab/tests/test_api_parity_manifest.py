# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Keep MATLAB-to-Python API classifications complete and importable."""

import sys
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS))

from audit_api_parity import audit


def test_matlab_api_manifest_is_complete():
    assert audit() == []
