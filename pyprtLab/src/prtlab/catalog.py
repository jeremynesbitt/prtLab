# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Reference optical systems used by examples and documentation."""

from __future__ import annotations

import numpy as np

from .systems import add_surface, create_optical_system


def example_cell_phone_lens_ex9_system():
    """Return the US7535658/Chipman cell-phone lens example 9."""

    radius = np.array(
        [1.745, -53.775, -2.133, 252.895, 4.001,
         165.253, 1.331, 1.269, np.inf, np.inf]
    )
    n_after = np.array(
        [1.471, 1.0, 1.606, 1.0, 1.510, 1.0, 1.510, 1.0, 1.516, 1.0]
    )
    asphere_data = np.array(
        [
            [1.955, -0.005327, -0.004451, -0.07406, 0.04328,
             0.02304, -0.02603, -0.06696, 0.0338],
            [9.972, -0.01241, -0.005103, -0.02493, -0.02322,
             0.0151, 0.0137, -0.05269, 0.003281],
            [1.962, -0.06145, -0.04302, 0.1481, -0.06121,
             -0.09732, 0.04933, 0.1549, -0.1473],
            [-0.0005579, -0.1786, -0.1257, 0.2152, 0.01437,
             -0.05593, -0.01349, 0.01755, 0.0007869],
            [-9.816, -0.1353, 0.09098, 0.004898, -0.01887,
             -0.01729, 0.003071, 0.0094, -0.005157],
            [-9.875, 0.04321, -0.01702, -0.01088, 7.168e-06,
             -0.001519, -0.0006694, 0.0009525, -0.000319],
            [-4.732, 0.04365, -0.09128, -0.003675, 0.006923,
             0.005455, -0.001346, -0.0003263, 6.607e-05],
            [-3.418, 0.04242, -0.1011, 0.04117, -0.006172,
             -0.002459, 0.0009779, 0.0001858, -9.145e-05],
            np.zeros(9),
            np.zeros(9),
        ]
    )
    thickness_after = np.array(
        [0.89, 0.69, 0.54, 0.07, 0.71, 0.09, 0.85, 0.70, 0.30,
         0.3799696104]
    )
    return _cell_phone_lens_system(
        "Chipman cell phone lens example 9",
        radius, n_after, asphere_data, thickness_after,
    )


def example_cell_phone_lens_ex3_system():
    """Return the US7535658/Chipman cell-phone lens example 3."""

    radius = np.array(
        [1.948, 9.999, -1.925, -6.492, 6.217,
         20.578, 1.399, 1.441, np.inf, np.inf]
    )
    n_after = np.array(
        [1.592, 1.0, 1.604, 1.0, 1.510, 1.0, 1.510, 1.0, 1.516, 1.0]
    )
    asphere_data = np.array(
        [
            [2.153, -3.226e-2, 7.995e-2, -9.382e-2, -4.859e-2,
             -1.676e-2, 8.074e-2, 1.276e-1, -1.544e-1],
            [40.18, -9.097e-3, -1.949e-2, -7.959e-3, 3.896e-4,
             1.469e-2, -2.080e-2, -9.120e-2, 6.076e-2],
            [2.105, -7.770e-2, -9.859e-3, 1.389e-1, -8.564e-2,
             -9.350e-2, 7.209e-2, 1.642e-1, -1.737e-1],
            [3.382, -1.389e-1, -1.206e-1, 2.242e-1, 1.777e-2,
             -6.023e-2, -1.833e-2, 1.595e-2, 1.589e-3],
            [-221.1, -1.030e-1, 1.165e-1, 1.276e-2, -2.206e-2,
             -1.973e-2, 2.674e-3, 1.031e-2, -3.570e-3],
            [0.9331, -1.121e-2, -1.583e-3, -4.034e-3, 8.135e-4,
             -1.640e-3, -4.422e-4, 1.210e-3, -3.105e-4],
            [-7.617, 1.131e-1, -1.116e-1, -1.101e-2, 5.557e-3,
             5.804e-3, -1.076e-3, -2.232e-4, 4.713e-5],
            [-2.707, 2.208e-2, -8.812e-2, 3.345e-2, -7.097e-3,
             -1.661e-3, 1.287e-3, 2.108e-4, -1.533e-4],
            np.zeros(9),
            np.zeros(9),
        ]
    )
    thickness_after = np.array(
        [0.9, 0.8, 0.59, 0.11, 0.66, 0.15, 0.7, 0.7, 0.3, 0.73]
    )
    return _cell_phone_lens_system(
        "Chipman cell phone lens example 3",
        radius, n_after, asphere_data, thickness_after,
    )


def _cell_phone_lens_system(
    description, radius, n_after, asphere_data, thickness_after
):
    system = create_optical_system(0.5876e-3)
    system.wavelength_units = "mm"
    system.description = description
    thickness_to_surface = np.r_[0, thickness_after[:-1]]
    n_before = np.r_[1.0, n_after[:-1]]
    bare = {"type": "bare"}
    for index in range(radius.size):
        surface_type = "plane" if np.isinf(radius[index]) else "odd_asphere"
        surface_data = {
            "kappa": asphere_data[index, 0],
            "A": asphere_data[index, 1:],
            "clearAperture": 1.0,
        }
        add_surface(
            system, radius[index], thickness_to_surface[index],
            surface_type, surface_data, "isotropic",
            {"n": n_before[index]}, {}, bare,
        )
    add_surface(
        system, np.inf, thickness_after[-1], "plane",
        {"note": "image plane"}, "isotropic", {"n": n_after[-1]}, {}, bare,
    )
    add_surface(
        system, np.inf, 0, "plane",
        {"note": "output medium after image plane"}, "isotropic",
        {"n": n_after[-1]}, {}, bare,
    )
    return system


def example_plane_mirror_system():
    """Return the minimal air-to-aluminum plane-mirror system."""

    system = create_optical_system(0.500e-3)
    system.wavelength_units = "mm"
    system.description = "Plane mirror"
    bare = {"type": "bare"}
    add_surface(
        system, np.inf, 0, "plane",
        {"note": "plane mirror", "clearAperture": 1.0},
        "isotropic", {"n": 1.0}, {}, bare,
    )
    add_surface(
        system, np.inf, 0, "plane",
        {"note": "complex-index metal side", "clearAperture": 1.0},
        "isotropic", {"n": 0.958 + 6.69j}, {}, bare,
    )
    return system


def example_cassegrain_telescope():
    """Return the two-mirror Cassegrain system from Chipman Table 12.7."""
    wavelength = 0.500e-3
    system = create_optical_system(wavelength)
    system.wavelength_units = "mm"
    system.description = "Cassegrain Telescope"
    mirror_index = 0.958 + 6.69j
    radii = [-5688.0, -1228.28]
    kappas = [0.0, -0.78]
    thicknesses = [0.0, -2317.80]
    indices_before = [1.0, mirror_index]
    bare = {"type": "bare"}
    for radius, kappa, thickness, index_before in zip(
        radii, kappas, thicknesses, indices_before, strict=True
    ):
        add_surface(
            system, radius, thickness, "odd_asphere",
            {"kappa": kappa, "A": np.zeros(8), "clearAperture": 4161.65},
            "isotropic", {"n": index_before}, {}, bare,
        )
    add_surface(
        system, np.inf, 3674.784, "plane", {"note": "image plane"},
        "isotropic", {"n": mirror_index}, {}, bare,
    )
    add_surface(
        system, np.inf, 0, "medium", {"note": "output medium"},
        "isotropic", {"n": 1.0}, {}, bare,
    )
    return system
