# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Reference optical-constant models used by prtLab examples."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def _return_like_input(value, reference):
    return float(value) if np.asarray(reference).ndim == 0 else value


def quartz_optical_constants(wavelength_um: ArrayLike):
    """Return ordinary and extraordinary quartz indices (wavelength in um)."""
    wavelength = np.asarray(wavelength_um, dtype=float)
    squared = wavelength**2
    n_o_squared = (
        1 + 0.28604141
        + 1.07044083 * squared / (squared - 1.00585997e-2)
        + 1.10202242 * squared / (squared - 100)
    )
    n_e_squared = (
        1 + 0.28851804
        + 1.09509924 * squared / (squared - 1.02101864e-2)
        + 1.15662475 * squared / (squared - 100)
    )
    return (
        _return_like_input(np.sqrt(n_o_squared), wavelength_um),
        _return_like_input(np.sqrt(n_e_squared), wavelength_um),
    )


def mgf2_optical_constants(wavelength_um: ArrayLike):
    """Return ordinary and extraordinary MgF2 indices (wavelength in um)."""
    wavelength = np.asarray(wavelength_um, dtype=float)
    squared = wavelength**2
    n_o_squared = (
        1 + 0.48755108 * squared / (squared - 0.04338408**2)
        + 0.39875031 * squared / (squared - 0.09461442**2)
        + 2.3120353 * squared / (squared - 23.793604**2)
    )
    n_e_squared = (
        1 + 0.41344023 * squared / (squared - 0.03684262**2)
        + 0.50497499 * squared / (squared - 0.09076162**2)
        + 2.4904862 * squared / (squared - 23.771995**2)
    )
    return (
        _return_like_input(np.sqrt(n_o_squared), wavelength_um),
        _return_like_input(np.sqrt(n_e_squared), wavelength_um),
    )


def n_pk52a_optical_constants(wavelength_um: ArrayLike):
    """Return the isotropic Schott N-PK52A index (wavelength in um)."""
    wavelength = np.asarray(wavelength_um, dtype=float)
    squared = wavelength**2
    n_squared = (
        1 + 1.029607 * squared / (squared - 0.00516800155)
        + 0.1880506 * squared / (squared - 0.0166658798)
        + 0.736488165 * squared / (squared - 138.964129)
    )
    return _return_like_input(np.sqrt(n_squared), wavelength_um)
