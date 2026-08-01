# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""MATLAB/Python parity for 3D P-to-Jones coordinate transforms."""

from __future__ import annotations

import numpy as np
import pytest
from parity import load_json_baseline

from prtlab.analysis import coherent_final_jones
from prtlab.coordinates import (
    _jones_local_basis,
    transform_p_to_jones,
)
from prtlab.models import RayBranch, TraceOptions, TraceResult

BASELINE = load_json_baseline("coordinateTransformBaselines.json")
TOLERANCE = BASELINE["absoluteTolerance"]


def unpack(value):
    return np.asarray(value["real"]) + 1j * np.asarray(value["imag"])


def coordinates(value):
    return {
        key: item if key == "type" else unpack(item)
        for key, item in value.items()
    }


@pytest.mark.parametrize(
    "reference", BASELINE["cases"], ids=lambda value: value["name"]
)
def test_transform_p_to_jones_matches_matlab(reference):
    p_matrix = unpack(reference["P"])
    k_in = unpack(reference["ki"])
    k_out = unpack(reference["kj"])
    coord_in = coordinates(reference["coordIn"])
    coord_out = coordinates(reference["coordOut"])
    actual = transform_p_to_jones(
        p_matrix, k_in, k_out, coord_in, coord_out
    )
    np.testing.assert_allclose(actual, unpack(reference["J"]), atol=TOLERANCE)
    x_out, y_out = _jones_local_basis(k_out, coord_out)
    np.testing.assert_allclose(x_out, unpack(reference["xOut"]), atol=TOLERANCE)
    np.testing.assert_allclose(y_out, unpack(reference["yOut"]), atol=TOLERANCE)


def test_output_coordinates_default_to_input_coordinates():
    reference = BASELINE["cases"][0]
    actual = transform_p_to_jones(
        unpack(reference["P"]), unpack(reference["ki"]),
        unpack(reference["kj"]), coordinates(reference["coordIn"]),
    )
    np.testing.assert_allclose(actual, unpack(reference["J"]), atol=TOLERANCE)


def test_sp_coordinates_reject_normal_incidence_without_a_basis_convention():
    with pytest.raises(ValueError, match="singular"):
        transform_p_to_jones(
            np.eye(3), [0, 0, 1], [0, 0, 1],
            {"type": "sp", "normal": [0, 0, 1]},
        )


def test_unknown_coordinate_type_is_rejected():
    with pytest.raises(ValueError, match="unknown coordinate type"):
        transform_p_to_jones(
            np.eye(3), [0, 0, 1], [0, 0, 1], {"type": "mystery"}
        )


def test_coherent_final_jones_accepts_distinct_output_coordinates():
    reference = BASELINE["cases"][1]
    ray = RayBranch(id=1, p_matrix=unpack(reference["P"]))
    result = TraceResult(None, TraceOptions(), rays=[ray], final_ray_ids=[1])
    actual = coherent_final_jones(
        result, unpack(reference["ki"]), unpack(reference["kj"]),
        coordinates(reference["coordIn"]), coordinates(reference["coordOut"]),
    )
    np.testing.assert_allclose(actual, unpack(reference["J"]), atol=TOLERANCE)


def test_missing_coordinate_field_has_actionable_error():
    with pytest.raises(ValueError, match="requires field 'a_loc'"):
        transform_p_to_jones(
            np.eye(3), [0, 0, 1], [0, 0, 1], {"type": "dipole"}
        )
