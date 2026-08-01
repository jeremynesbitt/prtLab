# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Trace an oblique polarized ray from air onto bare aluminum."""

import numpy as np

from prtlab import example_plane_mirror_system, polarization_ray_trace


def run_example(angle_degrees=30.0):
    angle = np.deg2rad(angle_degrees)
    system = example_plane_mirror_system()
    result = polarization_ray_trace(
        system,
        [0, np.sin(angle), np.cos(angle)],
        [0, 0, 0],
        np.array([1, 1]) / np.sqrt(2),
    )
    interaction = result.interactions[0]
    reflected = next(
        child for child in interaction.children
        if child.branch_type == "reflected"
    )
    return system, result, reflected


if __name__ == "__main__":
    _, _, ray = run_example()
    print("Reflected flux:", ray.flux)
    print("Reflected field:", ray.field_e)
    print("Accumulated P matrix:\n", ray.P)
