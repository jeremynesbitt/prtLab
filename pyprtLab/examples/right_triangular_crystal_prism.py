# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Construct and display a calcite right-triangular prism."""

import matplotlib.pyplot as plt

from prtlab import create_right_triangular_prism_solid, plot_prt_system_3d


def build_example():
    return create_right_triangular_prism_solid(
        {"nO": 1.6557, "nE": 1.4852},
        leg_y=12, leg_z=12, width=8, optic_axis=[1, 0, 0],
        wavelength=0.633, wavelength_units="um",
        name="Calcite right triangular prism",
    )


if __name__ == "__main__":
    prism = build_example()
    plot_prt_system_3d(
        prism, None, {"ShowRays": False, "ShowPolarization": False}
    )
    plt.show()
