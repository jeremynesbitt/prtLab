# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Quarter-wave-plate system construction for the Python port."""

from __future__ import annotations

import numpy as np

from prtlab import (
    add_surface,
    create_optical_system,
    polarization_ray_trace,
)


def example_quarter_wave_plate_system():
    wavelength = 0.633
    n_o = 1.656
    n_e = 1.485
    thickness = 1.25 * wavelength / (n_o - n_e)
    optic_axis = np.array([-np.sqrt(0.5), -np.sqrt(0.5), 0.0])

    system = create_optical_system(wavelength)
    system.wavelength_units = "um"
    bare = {"type": "bare"}
    add_surface(
        system, np.inf, 0, "plane", {}, "isotropic",
        {"n": 1.0}, {}, bare,
    )
    add_surface(
        system, np.inf, thickness, "plane", {}, "uniaxial",
        {"nO": n_o, "nE": n_e}, {"opticAxis": optic_axis}, bare,
    )
    add_surface(
        system, np.inf, 0, "plane", {}, "isotropic",
        {"n": 1.0}, {}, bare,
    )
    return system


if __name__ == "__main__":
    output = polarization_ray_trace(
        example_quarter_wave_plate_system(),
        [0, 0, 1],
        [0, 0, 0],
        [0, 1],
        {"encodePropagationPhaseInP": True},
    )
    p_total = sum(
        (output.rays[ray_id - 1].P for ray_id in output.final_ray_ids),
        np.zeros((3, 3), dtype=complex),
    )
    eigenvalues = np.linalg.eigvals(p_total[:2, :2])
    retardance = np.mod(
        np.angle(eigenvalues[1] / eigenvalues[0]) / (2 * np.pi), 1
    )
    print("Final forward ray IDs:", output.final_ray_ids)
    print("Retardance [waves]:", min(retardance, 1 - retardance))
