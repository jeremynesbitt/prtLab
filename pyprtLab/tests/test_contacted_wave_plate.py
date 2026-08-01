# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import sys
from pathlib import Path

import numpy as np

from prtlab import polarization_ray_trace

EXAMPLES = Path(__file__).parents[1] / "examples"
sys.path.insert(0, str(EXAMPLES))
from contacted_wave_plate import example_contacted_wave_plate_system


def test_contacted_plate_matches_matlab_accumulated_p():
    output = polarization_ray_trace(
        example_contacted_wave_plate_system(),
        [0, 0, 1],
        [0, 0, 0],
        np.array([1, 1]) / np.sqrt(2),
        {"encodePropagationPhaseInP": True},
    )
    assert len(output.final_ray_ids) == 2
    p_total = sum(
        (output.rays[ray_id - 1].P for ray_id in output.final_ray_ids),
        np.zeros((3, 3), dtype=complex),
    )
    expected = np.array(
        [
            [
                0.94046677985269 + 0.126783807935917j,
                5.75869815825507e-17 + 7.7632692286217e-18j,
                0,
            ],
            [
                -4.32582194212173e-19 + 5.81062962986831e-17j,
                0.0070646033536094 - 0.948947832781484j,
                0,
            ],
            [
                -3.18580156421384e-18 - 8.53264524022772e-18j,
                0,
                2,
            ],
        ],
        dtype=complex,
    )
    np.testing.assert_allclose(p_total, expected, atol=2e-12)
    np.testing.assert_allclose(
        [output.rays[ray_id - 1].opl for ray_id in output.final_ray_ids],
        [9.3375, 9.5085],
        atol=1e-12,
    )
    np.testing.assert_allclose(
        [output.rays[ray_id - 1].flux for ray_id in output.final_ray_ids],
        [0.450275948980610, 0.450275948980609],
        atol=1e-12,
    )


def test_uniaxial_interface_boundary_conditions():
    output = polarization_ray_trace(
        example_contacted_wave_plate_system(),
        [0, 0, 1],
        [0, 0, 0],
        np.array([1, 1]) / np.sqrt(2),
    )
    interactions = [
        item
        for item in output.interactions
        if item.case_name == "uniaxialToUniaxial"
    ]
    assert len(interactions) == 2
    for interaction in interactions:
        np.testing.assert_allclose(
            interaction.diagnostics["boundaryResidual_m"], 0, atol=1e-12
        )
