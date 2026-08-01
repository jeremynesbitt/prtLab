# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import numpy as np
from parity import load_json_baseline

from prtlab import example_cassegrain_telescope, polarization_ray_trace


def test_cassegrain_marginal_ray_matches_matlab():
    baseline = load_json_baseline("cassegrainMarginalRayBaseline.json")
    input_data = baseline["input"]
    result = polarization_ray_trace(
        example_cassegrain_telescope(),
        input_data["k"], input_data["position"], input_data["field"],
        {"minAmplitude": input_data["minAmplitude"]},
    )
    assert len(result.rays) == 4
    assert len(result.interactions) == 3
    assert result.final_ray_ids == [4]
    positions = np.array(
        [interaction.position.real for interaction in result.interactions]
    )
    np.testing.assert_allclose(positions, baseline["interactionPositions"], atol=2e-10)
    final = result.rays[3]
    expected = baseline["final"]
    np.testing.assert_allclose(final.k.real, expected["k"], atol=2e-14)
    np.testing.assert_allclose(final.flux, expected["flux"], atol=2e-14)
    np.testing.assert_allclose(final.opl, expected["opl"], atol=2e-12)
    expected_p = np.array(expected["pReal"]) + 1j * np.array(expected["pImag"])
    np.testing.assert_allclose(final.P, expected_p, atol=4e-15)


def test_positive_thickness_system_does_not_advance_backward_reflection():
    from prtlab import add_surface, create_optical_system

    system = create_optical_system(0.55)
    add_surface(system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.0}, {})
    add_surface(system, np.inf, 1, "plane", {}, "isotropic", {"n": 1.5}, {})
    add_surface(system, np.inf, 1, "plane", {}, "isotropic", {"n": 1.0}, {})
    result = polarization_ray_trace(system, [0, 0, 1], [0, 0, 0], [1, 0])
    assert len(result.interactions) == 2
    assert len(result.final_ray_ids) == 1
    assert result.rays[result.final_ray_ids[0] - 1].branch_type == "transmitted"
