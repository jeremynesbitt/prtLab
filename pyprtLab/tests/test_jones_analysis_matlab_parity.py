# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""MATLAB/Python parity for Jones polarization ellipse analysis."""

import numpy as np
import pytest
from parity import load_json_baseline

from prtlab.analysis import jones_major_axis_ellipticity

BASELINE = load_json_baseline("jonesAnalysisBaselines.json")


def unpack(value):
    return np.asarray(value["real"]) + 1j * np.asarray(value["imag"])


def test_jones_ellipse_parameters_match_matlab():
    actual = jones_major_axis_ellipticity(unpack(BASELINE["fields"]))
    keys = ("psi", "ellipticity", "majorAxis", "major", "minor")
    for value, key in zip(actual, keys, strict=True):
        np.testing.assert_allclose(
            value, unpack(BASELINE[key]), atol=BASELINE["absoluteTolerance"]
        )


def test_single_circular_state_has_unit_ellipticity():
    psi, ellipticity, _, major, minor = jones_major_axis_ellipticity(
        np.array([1, 1j]) / np.sqrt(2)
    )
    np.testing.assert_allclose(ellipticity, 1, atol=1e-14)
    np.testing.assert_allclose(major, minor, atol=1e-14)
    assert np.isfinite(psi)


def test_jones_ellipse_analysis_rejects_wrong_shape():
    with pytest.raises(ValueError, match="shape"):
        jones_major_axis_ellipticity(np.ones((3, 2)))
