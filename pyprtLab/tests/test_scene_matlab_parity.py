# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""MATLAB/Python parity for branched solid and multi-solid traces."""

from __future__ import annotations

import numpy as np
import pytest
from parity import load_json_baseline

from prtlab import (
    create_fresnel_rhomb_solid,
    create_glan_taylor_polarizer_scene,
    jones_retardance,
    polarization_ray_trace_scene,
    polarization_ray_trace_solid,
)
from prtlab.coordinates import transform_p_to_jones

BASELINE = load_json_baseline("sceneTraceBaselines.json")
TOLERANCE = BASELINE["absoluteTolerance"]


def unpack(value):
    return np.asarray(value["real"]) + 1j * np.asarray(value["imag"])


def string_list(value):
    return [value] if isinstance(value, str) else value


def ray_path(result, ray):
    by_parent = {item.incident.id: item for item in result.interactions}
    interactions = [
        by_parent[ray_id]
        for ray_id in ray.metadata["history"]
        if ray_id in by_parent
    ]
    intercepts = [item.diagnostics["interceptData"] for item in interactions]
    return {
        "faces": [item["faceName"] for item in intercepts],
        "directions": [item["direction"] for item in intercepts],
        "interfaceCases": [item.case_name for item in interactions],
        "incidentModes": [item.incident.mode for item in interactions],
    }


@pytest.mark.parametrize(
    "reference", BASELINE["fresnelRhomb"],
    ids=lambda value: f"aoi-{value['angleDeg']}-deg",
)
def test_fresnel_rhomb_dominant_path_matches_matlab(reference):
    solid = create_fresnel_rhomb_solid(
        1.4965, length=14, height=8, width=8, shear=11.02085,
        wavelength=0.589,
    )
    angle = np.deg2rad(reference["angleDeg"])
    k_in = np.array([0, np.sin(angle), np.cos(angle)])
    result = polarization_ray_trace_solid(
        solid, k_in, [0, 5, -2], np.array([1, 1]) / np.sqrt(2),
        {"maxInteractions": 8, "minAmplitude": 1e-4},
    )
    ray = max(
        (result.rays[ray_id - 1] for ray_id in result.final_ray_ids),
        key=lambda item: abs(item.flux),
    )

    assert len(result.rays) == reference["rayCount"]
    assert len(result.interactions) == reference["interactionCount"]
    assert len(result.final_ray_ids) == reference["finalRayCount"]
    assert ray.id == reference["dominantRayId"]
    assert_ray_matches(result, ray, reference["dominant"])

    coordinates = {
        "type": "doublePole", "a_loc": np.array([0, 0, 1]),
        "x_o": np.array([1, 0, 0]),
    }
    jones = transform_p_to_jones(ray.P, k_in, ray.k, coordinates)
    assert_complex_close(jones, reference["jones"])
    np.testing.assert_allclose(
        jones_retardance(jones), unpack(reference["retardance"]),
        atol=TOLERANCE,
    )


def test_glan_taylor_surviving_paths_match_matlab():
    leg_y = 10.0
    scene = create_glan_taylor_polarizer_scene(
        {"nO": 1.656, "nE": 1.485}, leg_y=leg_y,
        leg_z=leg_y * np.tan(np.deg2rad(40)), width=8, air_gap=0.2,
        optic_axis=[1, 0, 0], wavelength=0.633,
        name="Calcite Glan-Taylor polarizer",
    )
    result = polarization_ray_trace_scene(
        scene, [0, 0, 1], [0, leg_y / 2, -1],
        np.array([1, 1]) / np.sqrt(2),
        {"maxInteractions": 8, "minAmplitude": 0.2,
         "minRelativeFlux": 1e-3},
    )
    reference = BASELINE["glanTaylor"]
    assert len(result.rays) == reference["rayCount"]
    assert len(result.interactions) == reference["interactionCount"]
    assert result.final_ray_ids == reference["finalRayIds"]
    for expected in reference["finalRays"]:
        assert_ray_matches(result, result.rays[expected["id"] - 1], expected)


def assert_ray_matches(result, ray, reference):
    assert ray.id == reference["id"]
    assert ray.mode == reference["mode"]
    assert ray.branch_type == reference["branchType"]
    assert ray.metadata["history"] == reference["history"]
    path = ray_path(result, ray)
    for key in ("faces", "directions", "interfaceCases", "incidentModes"):
        assert path[key] == string_list(reference[key])
    assert_complex_close(ray.position, reference["position"])
    assert_complex_close(ray.k, reference["k"])
    assert_complex_close(ray.S, reference["S"])
    assert_complex_close(ray.field_e, reference["fieldE"])
    assert_complex_close(ray.field_h, reference["fieldH"])
    assert_complex_close(ray.P, reference["P"])
    assert_complex_close(ray.Q, reference["Q"])
    np.testing.assert_allclose(ray.amplitude, unpack(reference["amplitude"]), atol=TOLERANCE)
    np.testing.assert_allclose(ray.flux, unpack(reference["flux"]), atol=TOLERANCE)
    np.testing.assert_allclose(ray.opl, unpack(reference["OPL"]), atol=TOLERANCE)


def assert_complex_close(actual, expected):
    np.testing.assert_allclose(actual, unpack(expected), atol=TOLERANCE)
