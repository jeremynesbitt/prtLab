# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import prtlab


def test_runtime_version_matches_distribution_metadata():
    from importlib.metadata import version

    assert prtlab.__version__ == version("prtLab")


def test_public_api_is_explicit_and_importable():
    assert len(prtlab.__all__) == len(set(prtlab.__all__))
    for name in prtlab.__all__:
        assert hasattr(prtlab, name), name
    assert prtlab.__version__ == "0.1.0.dev0"
