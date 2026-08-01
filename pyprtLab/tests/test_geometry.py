# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import numpy as np

from prtlab.geometry import intersect_surface, odd_asphere_sag
from prtlab.models import RayBranch, Surface


def test_spherical_conic_sag():
    radius = 100.0
    radial_coordinate = 2.0
    sag = odd_asphere_sag(
        radial_coordinate, 0, radius, 1.0, np.zeros(8)
    )
    expected = radius - np.sqrt(radius**2 - radial_coordinate**2)
    np.testing.assert_allclose(sag, expected, atol=1e-14)


def test_odd_polynomial_term():
    coefficients = np.zeros(8)
    coefficients[0] = 2e-4
    sag = odd_asphere_sag(3, 4, np.inf, 0, coefficients)
    np.testing.assert_allclose(sag, 2e-4 * 5**3)


def test_vertical_asphere_intercept_satisfies_sag_equation():
    surface = Surface(
        radius=50,
        thickness=10,
        surface_type="odd_asphere",
        surface_data={"kappa": 1.0, "A": [1e-6]},
        material_type="isotropic",
        index_data={"n": 1.0},
    )
    ray = RayBranch(
        position=np.array([0.0, 2.0, 0.0]),
        k=np.array([0.0, 0.0, 1.0]),
        s_direction=np.array([0.0, 0.0, 1.0]),
        metadata={"currentVertexZ": 0.0},
    )
    hit, normal, data = intersect_surface(surface, ray)
    expected_sag = odd_asphere_sag(
        hit[0].real, hit[1].real, surface.radius, 1.0,
        np.array([1e-6, 0, 0, 0, 0, 0, 0, 0]),
    )
    np.testing.assert_allclose(hit[2].real, 10 + expected_sag, atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(normal), 1, atol=1e-13)
    assert abs(data["residual"]) < 1e-11
