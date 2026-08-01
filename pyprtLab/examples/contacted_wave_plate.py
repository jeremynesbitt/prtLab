# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Two optically contacted uniaxial plates."""

from __future__ import annotations

import numpy as np

from prtlab import add_surface, create_optical_system


def example_contacted_wave_plate_system():
    wavelength = 0.633
    n_o = 1.656
    n_e = 1.485
    system = create_optical_system(wavelength)
    system.wavelength_units = "um"
    bare = {"type": "bare"}
    add_surface(
        system, np.inf, 0, "plane", {}, "isotropic",
        {"n": 1.0}, {}, bare,
    )
    add_surface(
        system, np.inf, 2.5, "plane", {}, "uniaxial",
        {"nO": n_o, "nE": n_e}, {"opticAxis": np.array([1.0, 0.0, 0.0])},
        bare,
    )
    add_surface(
        system, np.inf, 3.5, "plane", {}, "uniaxial",
        {"nO": n_o, "nE": n_e}, {"opticAxis": np.array([0.0, 1.0, 0.0])},
        bare,
    )
    add_surface(
        system, np.inf, 0, "plane", {}, "isotropic",
        {"n": 1.0}, {}, bare,
    )
    return system
