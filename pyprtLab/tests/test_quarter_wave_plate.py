# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import sys
from pathlib import Path

import numpy as np

from prtlab import polarization_ray_trace

EXAMPLES = Path(__file__).parents[1] / "examples"
sys.path.insert(0, str(EXAMPLES))
from quarter_wave_plate import example_quarter_wave_plate_system


def test_quarter_wave_plate_matches_matlab_accumulated_p():
    output = polarization_ray_trace(
        example_quarter_wave_plate_system(),
        [0, 0, 1],
        [0, 0, 0],
        [0, 1],
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
                0.6659084104843634 - 0.09116845886285108j,
                -0.07509212227185058 - 0.6679123237506079j,
                0,
            ],
            [
                -0.07509212227185069 - 0.6679123237506079j,
                0.6659084104843632 - 0.09116845886285158j,
                0,
            ],
            [0, 0, 2],
        ],
        dtype=complex,
    )
    np.testing.assert_allclose(p_total, expected, atol=2e-12)
    np.testing.assert_allclose(
        [output.rays[ray_id - 1].opl for ray_id in output.final_ray_ids],
        [7.662631578947377, 6.871381578947378],
        atol=1e-12,
    )

    eigenvalues = np.linalg.eigvals(p_total[:2, :2])
    retardance = np.mod(
        np.angle(eigenvalues[1] / eigenvalues[0]) / (2 * np.pi), 1
    )
    np.testing.assert_allclose(min(retardance, 1 - retardance), 0.25,
                               atol=1e-12)


def test_propagation_phase_is_disabled_by_default():
    output = polarization_ray_trace(
        example_quarter_wave_plate_system(),
        [0, 0, 1],
        [0, 0, 0],
        [0, 1],
    )
    assert not output.options.encode_propagation_phase_in_p
    p_total = sum(
        (output.rays[ray_id - 1].P for ray_id in output.final_ray_ids),
        np.zeros((3, 3), dtype=complex),
    )
    np.testing.assert_allclose(p_total.imag, 0, atol=1e-13)
    assert all(
        ray.metadata.get("lastSegmentPhase") is None for ray in output.rays
    )
