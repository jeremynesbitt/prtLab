# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import numpy as np

from prtlab import add_surface, create_optical_system, polarization_ray_trace
from prtlab.linalg import calc_k_from_theta_x_theta_y


def glass_interface():
    system = create_optical_system(0.55)
    add_surface(
        system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.0}, {}
    )
    add_surface(
        system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.5}, {}
    )
    return system


def test_normal_incidence_coefficients_match_fresnel():
    result = polarization_ray_trace(
        glass_interface(), [0, 0, 1], [0, 0, 0], [1, 0]
    )
    interaction = result.interactions[0]
    expected_t = 2 / 2.5
    expected_r = (1 - 1.5) / 2.5
    np.testing.assert_allclose(interaction.coefficients["As"][0], expected_t)
    np.testing.assert_allclose(
        abs(interaction.coefficients["As"][2]), abs(expected_r)
    )
    np.testing.assert_allclose(
        interaction.diagnostics["boundaryResidual_s"], 0, atol=1e-13
    )
    np.testing.assert_allclose(
        interaction.diagnostics["boundaryResidual_p"], 0, atol=1e-13
    )


def test_oblique_interface_conserves_tangential_q():
    k_in = calc_k_from_theta_x_theta_y(20, 15)
    result = polarization_ray_trace(
        glass_interface(), k_in, [0, 0, 0], [1, 0]
    )
    transmitted = next(
        child
        for child in result.interactions[0].children
        if child.branch_type == "transmitted"
    )
    np.testing.assert_allclose(
        1.5 * transmitted.k[:2], k_in[:2], atol=1e-13
    )
