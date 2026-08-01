# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Trace ordinary and extraordinary paths through a Glan-Taylor polarizer."""

import numpy as np

from prtlab import (
    create_glan_taylor_polarizer_scene,
    plot_prt_system_3d,
    polarization_ray_trace_scene,
)


def run_example():
    leg_y = 10.0
    scene = create_glan_taylor_polarizer_scene(
        {"nO": 1.656, "nE": 1.485}, leg_y=leg_y,
        leg_z=leg_y * np.tan(np.deg2rad(40)), width=8, air_gap=0.2,
        optic_axis=[1, 0, 0], wavelength=0.633,
        name="Calcite Glan-Taylor polarizer",
    )
    output = polarization_ray_trace_scene(
        scene, [0, 0, 1], [0, leg_y / 2, -1], [1, 1] / np.sqrt(2),
        {"maxInteractions": 8, "minAmplitude": 0.2, "minRelativeFlux": 1e-3},
    )
    return scene, output


def plot_example():
    scene, output = run_example()
    return plot_prt_system_3d(
        scene, output,
        {"PreExtend": 1, "PostExtend": 1, "AxisPaddingFraction": 0.2},
    )


if __name__ == "__main__":
    _, result = run_example()
    for ray_id in result.final_ray_ids:
        ray = result.rays[ray_id - 1]
        print(
            f"ray {ray_id}: {ray.mode:13s} {ray.branch_type:11s} "
            f"flux={ray.flux:.8g}"
        )
