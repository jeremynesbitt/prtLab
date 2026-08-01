# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Detailed MATLAB/Python parity checks for anisotropic interfaces."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from parity import load_json_baseline

from prtlab import add_surface, create_optical_system
from prtlab.interfaces import (
    trace_isotropic_to_uniaxial,
    trace_uniaxial_to_isotropic,
    trace_uniaxial_to_uniaxial,
)
from prtlab.tracing import _initial_ray

BASELINE_NAME = "anisotropicInterfaceBaselines.json"
BASELINE_PATH = Path(__file__).parent / "baselines" / BASELINE_NAME
BASELINE = load_json_baseline(BASELINE_NAME)


def unpack(value):
    """Decode the explicit complex-number representation emitted by MATLAB."""

    return np.asarray(value["real"]) + 1j * np.asarray(value["imag"])


def make_system(case, final_material):
    system = create_optical_system(BASELINE["wavelength"])
    bare = {"type": "bare"}
    add_surface(
        system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.0}, {}, bare
    )
    add_surface(
        system, np.inf, 1, "plane", {}, "uniaxial", case["index1"],
        {"opticAxis": unpack(case["axis1"]).real}, bare,
    )
    if final_material == "isotropic":
        add_surface(
            system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.0}, {}, bare
        )
    else:
        add_surface(
            system, np.inf, 0, "plane", {}, "uniaxial", case["index2"],
            {"opticAxis": unpack(case["axis2"]).real}, bare,
        )
    return system


def interactions_for_case(case, final_material):
    system = make_system(case, final_material)
    initial = _initial_ray(
        system, unpack(case["kIncident"]).real, [0, 0, 0],
        unpack(case["inputField"]),
    )
    entry = trace_isotropic_to_uniaxial(
        initial, system.surfaces[0], system.surfaces[1],
        np.zeros(3), np.array([0, 0, 1]), 1,
    )
    kernel = (
        trace_uniaxial_to_isotropic
        if final_material == "isotropic"
        else trace_uniaxial_to_uniaxial
    )
    return [
        kernel(
            incident, system.surfaces[1], system.surfaces[2],
            np.array([0, 0, 1]), np.array([0, 0, 1]), 2,
        )
        for incident in entry.children[:2]
    ]


@pytest.mark.parametrize("case", BASELINE["cases"], ids=lambda case: case["name"])
@pytest.mark.parametrize(
    ("fixture_key", "final_material"),
    [
        ("uniaxialToIsotropic", "isotropic"),
        ("uniaxialToUniaxial", "uniaxial"),
    ],
)
def test_anisotropic_interface_matches_matlab(case, fixture_key, final_material):
    tolerance = BASELINE["absoluteTolerance"]
    actual_interactions = interactions_for_case(case, final_material)
    expected_interactions = case[fixture_key]
    assert len(actual_interactions) == len(expected_interactions) == 2

    for actual, expected in zip(actual_interactions, expected_interactions):
        assert actual.incident.mode == expected["incidentMode"]
        assert_complex_close(actual.incident.k, expected["incidentK"], tolerance)
        assert_complex_close(actual.incident.s_direction, expected["incidentS"], tolerance)
        assert_complex_close(actual.incident.field_e, expected["incidentFieldE"], tolerance)
        assert_complex_close(actual.incident.field_h, expected["incidentFieldH"], tolerance)
        assert len(actual.children) == len(expected["children"])

        for child, reference in zip(actual.children, expected["children"]):
            assert child.mode == reference["mode"]
            assert child.branch_type == reference["branchType"]
            assert child.active == reference["active"]
            assert_complex_close(child.k, reference["k"], tolerance)
            assert_complex_close(child.s_direction, reference["S"], tolerance)
            assert_mode_pair_close(child, reference, tolerance)
            assert_complex_close(child.field_e, reference["fieldE"], tolerance)
            assert_complex_close(child.field_h, reference["fieldH"], tolerance)
            assert_complex_close(child.p_matrix, reference["P"], tolerance)
            assert_complex_close(child.metadata["P_interface"], reference["interfaceP"], tolerance)
            np.testing.assert_allclose(child.amplitude, unpack(reference["amplitude"]), atol=tolerance)
            np.testing.assert_allclose(child.flux, unpack(reference["flux"]), atol=tolerance)
            np.testing.assert_allclose(child.metadata["n"], unpack(reference["index"]), atol=tolerance)

        np.testing.assert_allclose(
            actual.diagnostics["boundaryResidual_m"], 0, atol=1e-12
        )
        np.testing.assert_allclose(
            sum(child.flux for child in actual.children),
            actual.incident.flux,
            atol=5e-10,
        )


def assert_mode_pair_close(child, reference, tolerance):
    expected_e = unpack(reference["modeE"])
    expected_h = unpack(reference["modeH"])
    overlap = np.vdot(expected_e, child.mode_e)
    phase = 1 if abs(overlap) == 0 else overlap / abs(overlap)
    np.testing.assert_allclose(child.mode_e, phase * expected_e, atol=tolerance)
    np.testing.assert_allclose(child.mode_h, phase * expected_h, atol=tolerance)


def assert_complex_close(actual, expected, tolerance):
    np.testing.assert_allclose(actual, unpack(expected), atol=tolerance)
