# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Optimize a quartz/MgF2 broadband quarter-wave retarder."""

import matplotlib.pyplot as plt
import numpy as np

from prtlab import (
    achromatic_retardance_spectrum,
    optimize_achromatic_wave_plate,
    scalar_crossed_axis_retardance,
)


def run_example():
    wavelengths = np.linspace(0.4, 0.7, 7)
    optimization = optimize_achromatic_wave_plate(wavelengths)
    retardance = achromatic_retardance_spectrum(
        wavelengths, *optimization.x
    )
    scalar = scalar_crossed_axis_retardance(wavelengths, *optimization.x)
    return wavelengths, optimization, retardance, scalar


def plot_example():
    wavelengths, optimization, retardance, scalar = run_example()
    _, axis = plt.subplots()
    axis.plot(1000 * wavelengths, retardance, "o-", label="prtLab")
    axis.plot(1000 * wavelengths, scalar, "s:", label="scalar estimate")
    axis.axhline(0.25, color="black", linestyle="--", label="target")
    axis.set(xlabel="Wavelength (nm)", ylabel="Retardance (waves)")
    axis.grid(True)
    axis.legend()
    axis.set_title(
        f"Quartz {optimization.x[0]:.3f} um / MgF2 {optimization.x[1]:.3f} um"
    )
    return axis


if __name__ == "__main__":
    _, result, retardance, _ = run_example()
    print("Quartz thickness (um):", result.x[0])
    print("MgF2 thickness (um):", result.x[1])
    print("Retardance (waves):", retardance)
