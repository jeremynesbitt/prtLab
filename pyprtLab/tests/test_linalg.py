# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import numpy as np

from prtlab.linalg import (
    calc_k_from_theta_x_theta_y,
    get_theta_x_theta_y_from_k,
    make_k,
    normalize,
    snell_vector,
)


def test_make_k_is_cross_product_operator():
    a = np.array([0.2, -0.3, 0.9])
    b = np.array([-0.4, 0.7, 0.1])
    np.testing.assert_allclose(make_k(a) @ b, np.cross(a, b))


def test_angle_convention_round_trip_components():
    k = calc_k_from_theta_x_theta_y(30, -20)
    theta_x, theta_y = get_theta_x_theta_y_from_k(k)
    np.testing.assert_allclose(np.rad2deg(theta_x), 30, atol=1e-12)
    np.testing.assert_allclose(np.rad2deg(theta_y), -20, atol=1e-12)


def test_vector_snell_preserves_tangential_wavevector():
    k_in = calc_k_from_theta_x_theta_y(25, 10)
    normal = np.array([0.0, 0.0, 1.0])
    k_out = snell_vector(k_in, normal, 1.0, 1.5)
    np.testing.assert_allclose(
        1.5 * (k_out - np.dot(k_out, normal) * normal),
        k_in - np.dot(k_in, normal) * normal,
        atol=1e-13,
    )
    np.testing.assert_allclose(np.dot(k_out, k_out), 1.0, atol=1e-13)


def test_normalize_real_vector():
    value = normalize([3, 0, 4])
    np.testing.assert_allclose(value, [0.6, 0, 0.8])
