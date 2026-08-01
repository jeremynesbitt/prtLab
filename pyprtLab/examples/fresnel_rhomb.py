# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Trace a diagonal input polarization through a Fresnel rhomb."""

import numpy as np

from prtlab import (
    create_fresnel_rhomb_solid,
    plot_prt_system_3d,
    polarization_ray_trace_solid,
)


def run_example():
    solid = create_fresnel_rhomb_solid(
        1.495, length=14, height=8, width=8, shear=11.02085,
        wavelength=0.589,
    )
    return polarization_ray_trace_solid(
        solid,
        [0, np.sin(np.deg2rad(0.1)), np.cos(np.deg2rad(0.1))],
        [0, 5, -2], [1, 1] / np.sqrt(2),
        {"maxInteractions": 8, "minAmplitude": 1e-4},
    )


def plot_example():
    solid = create_fresnel_rhomb_solid(
        1.495, length=14, height=8, width=8, shear=11.02085,
        wavelength=0.589,
    )
    result = run_example()
    return plot_prt_system_3d(
        solid, result,
        {"PostExtend": 4, "AxisPaddingFraction": 0.22, "ShowLegend": True},
    )


if __name__ == "__main__":
    output = run_example()
    for interaction in output.interactions:
        data = interaction.diagnostics["interceptData"]
        print(f"{data['faceName']}: {data['direction']}")
