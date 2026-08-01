# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import numpy as np

from prtlab import (
    create_fresnel_rhomb_solid,
    create_glan_taylor_polarizer_scene,
    create_right_triangular_prism_solid,
    polarization_ray_trace_scene,
    polarization_ray_trace_solid,
    transform_solid,
)


def _all_transmitted_final(result):
    for ray_id in result.final_ray_ids:
        ray = result.rays[ray_id - 1]
        history = ray.metadata["history"]
        branches = [result.rays[index - 1].branch_type for index in history]
        if all(branch == "transmitted" for branch in branches[1:]):
            return ray
    raise AssertionError("no all-transmitted final ray was found")


def test_two_prism_air_gap_and_dominant_path():
    gap = 0.2
    scene = create_glan_taylor_polarizer_scene(
        {"n": 1.2}, leg_y=10, leg_z=10, width=8,
        air_gap=gap, wavelength=0.633,
    )
    h1 = next(face for face in scene.solids[0].faces if face.name == "hypotenuse")
    h2 = next(face for face in scene.solids[1].faces if face.name == "hypotenuse")
    np.testing.assert_allclose(abs(np.dot(h2.point - h1.point, h1.normal)), gap, atol=1e-12)

    result = polarization_ray_trace_scene(
        scene, [0, 0, 1], [0, 5, -1], [1, 0],
        {"maxInteractions": 8, "minAmplitude": 1e-8, "minRelativeFlux": 0.02},
    )
    ray = _all_transmitted_final(result)
    assert ray.metadata["solidIndex"] == 0
    history = ray.metadata["history"][:-1]
    interactions = [
        next(item for item in result.interactions if item.incident.id == ray_id)
        for ray_id in history
    ]
    intercepts = [item.diagnostics["interceptData"] for item in interactions]
    assert [item["faceName"] for item in intercepts] == [
        "z leg", "hypotenuse", "hypotenuse", "z leg"
    ]
    assert [item["solidIndex"] for item in intercepts] == [1, 1, 2, 2]
    assert [item["direction"] for item in intercepts] == [
        "entering", "exiting", "entering", "exiting"
    ]
    np.testing.assert_allclose(ray.position, [0, 4.914642672431773, 10.141421356237309], atol=1e-13)
    np.testing.assert_allclose(ray.k, [0, 0, 1], atol=2e-13)
    np.testing.assert_allclose(ray.flux, 0.8806652601406875, atol=1e-13)
    np.testing.assert_allclose(ray.opl, 13.105830052442585, atol=1e-13)
    np.testing.assert_allclose(
        ray.p_matrix,
        np.diag([0.9384376698218627, 0.9888712023351989, 1.0]),
        atol=2e-13,
    )


def test_single_solid_wrapper_traces_fresnel_rhomb():
    solid = create_fresnel_rhomb_solid(
        1.495, length=14, height=8, width=8, shear=11.02085,
        wavelength=0.589,
    )
    result = polarization_ray_trace_solid(
        solid,
        [0, np.sin(np.deg2rad(0.1)), np.cos(np.deg2rad(0.1))],
        [0, 5, -2], [1, 1] / np.sqrt(2),
        {"maxInteractions": 8, "minAmplitude": 1e-4},
    )
    assert result.interactions
    assert any(
        child.branch_type == "reflected"
        for interaction in result.interactions
        for child in interaction.children
    )


def test_rigid_transform_updates_faces_and_optic_axis():
    solid = create_right_triangular_prism_solid(
        {"nO": 1.656, "nE": 1.485}, leg_y=3, leg_z=4,
        width=2, optic_axis=[0, 1, 0],
    )
    rotation = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    translation = np.array([2, 3, 4])
    point = solid.faces[0].point.copy()
    normal = solid.faces[0].normal.copy()
    transformed = transform_solid(solid, rotation, translation)
    np.testing.assert_allclose(transformed.faces[0].point, rotation @ point + translation)
    np.testing.assert_allclose(transformed.faces[0].normal, rotation @ normal)
    np.testing.assert_allclose(transformed.material.axis_data["opticAxis"], [0, 0, 1])


def test_glan_taylor_branch_tree_matches_matlab():
    leg_y = 10.0
    scene = create_glan_taylor_polarizer_scene(
        {"nO": 1.656, "nE": 1.485}, leg_y=leg_y,
        leg_z=leg_y * np.tan(np.deg2rad(40)), width=8, air_gap=0.2,
        optic_axis=[1, 0, 0], wavelength=0.633,
    )
    result = polarization_ray_trace_scene(
        scene, [0, 0, 1], [0, leg_y / 2, -1], [1, 1] / np.sqrt(2),
        {"maxInteractions": 8, "minAmplitude": 0.2, "minRelativeFlux": 1e-3},
    )
    assert len(result.rays) == 15
    assert len(result.interactions) == 10
    assert result.final_ray_ids == [4, 8, 11, 12, 15]
    expected = [
        ("isotropic", "reflected", 0.049547388894541806),
        ("isotropic", "transmitted", 0.43707664914877492),
        ("isotropic", "transmitted", 0.15518932215423142),
        ("isotropic", "transmitted", 0.200349981298483),
        ("isotropic", "transmitted", 0.06720689333046544),
    ]
    for ray_id, (mode, branch, flux) in zip(result.final_ray_ids, expected, strict=True):
        ray = result.rays[ray_id - 1]
        assert (ray.mode, ray.branch_type) == (mode, branch)
        np.testing.assert_allclose(ray.flux, flux, atol=1e-13)
