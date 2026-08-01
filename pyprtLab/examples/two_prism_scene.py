# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Trace through two isotropic prisms separated by an air gap."""

from prtlab import create_glan_taylor_polarizer_scene, polarization_ray_trace_scene


def run_example():
    scene = create_glan_taylor_polarizer_scene(
        {"n": 1.2}, leg_y=10, leg_z=10, width=8,
        air_gap=0.2, wavelength=0.633, name="Two-prism scene test",
    )
    return polarization_ray_trace_scene(
        scene, [0, 0, 1], [0, 5, -1], [1, 0],
        {"maxInteractions": 8, "minAmplitude": 1e-8},
    )


if __name__ == "__main__":
    output = run_example()
    for interaction in output.interactions:
        data = interaction.diagnostics["interceptData"]
        print(f"{data['solidName']} / {data['faceName']} ({data['direction']})")
