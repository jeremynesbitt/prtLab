# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Constructors and broadband analysis for compound uniaxial waveplates."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .analysis import coherent_final_jones, jones_retardance_waves
from .dispersion import mgf2_optical_constants, quartz_optical_constants
from .systems import add_surface, create_optical_system
from .tracing import polarization_ray_trace

DEFAULT_QUARTZ_AXIS = np.array([-np.sqrt(0.5), -np.sqrt(0.5), 0.0])
DEFAULT_MGF2_AXIS = np.array([np.sqrt(0.5), -np.sqrt(0.5), 0.0])


def quarter_wave_plate_system(
    wavelength_um: float = 0.633,
    n_o: float = 1.656,
    n_e: float = 1.485,
    order: float = 1.25,
    optic_axis=DEFAULT_QUARTZ_AXIS,
):
    """Create the air/uniaxial/air plate used in Chipman Example 20.1."""
    thickness = order * wavelength_um / (n_o - n_e)
    system = create_optical_system(wavelength_um)
    system.wavelength_units = "um"
    bare = {"type": "bare"}
    add_surface(system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.0}, {}, bare)
    add_surface(
        system, np.inf, thickness, "plane", {}, "uniaxial",
        {"nO": n_o, "nE": n_e},
        {"opticAxis": np.asarray(optic_axis, dtype=float)}, bare,
    )
    add_surface(system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.0}, {}, bare)
    return system


def achromatic_wave_plate_system(
    wavelength_um: float,
    quartz_thickness: float = 230.537086,
    gap_thickness: float = 444.5,
    mgf2_thickness: float = 190.217463,
    quartz_axis=DEFAULT_QUARTZ_AXIS,
    mgf2_axis=DEFAULT_MGF2_AXIS,
):
    """Create an air/quartz/air/MgF2/air achromatic waveplate."""
    n_o_quartz, n_e_quartz = quartz_optical_constants(wavelength_um)
    n_o_mgf2, n_e_mgf2 = mgf2_optical_constants(wavelength_um)
    system = create_optical_system(wavelength_um)
    system.wavelength_units = "um"
    bare = {"type": "bare"}
    add_surface(system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.0}, {}, bare)
    add_surface(
        system, np.inf, quartz_thickness, "plane", {}, "uniaxial",
        {"nO": n_o_quartz, "nE": n_e_quartz},
        {"opticAxis": np.asarray(quartz_axis, dtype=float)}, bare,
    )
    add_surface(system, np.inf, gap_thickness, "plane", {}, "isotropic", {"n": 1.0}, {}, bare)
    add_surface(
        system, np.inf, mgf2_thickness, "plane", {}, "uniaxial",
        {"nO": n_o_mgf2, "nE": n_e_mgf2},
        {"opticAxis": np.asarray(mgf2_axis, dtype=float)}, bare,
    )
    add_surface(system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.0}, {}, bare)
    return system


def achromatic_retardance_spectrum(
    wavelengths_um,
    quartz_thickness: float,
    mgf2_thickness: float,
    *,
    gap_thickness: float = 444.5,
    quartz_axis=DEFAULT_QUARTZ_AXIS,
    mgf2_axis=DEFAULT_MGF2_AXIS,
):
    """Trace and return compound-waveplate retardance versus wavelength."""
    wavelengths = np.asarray(wavelengths_um, dtype=float)
    retardance = np.empty_like(wavelengths)
    for index, wavelength in np.ndenumerate(wavelengths):
        system = achromatic_wave_plate_system(
            float(wavelength), quartz_thickness, gap_thickness,
            mgf2_thickness, quartz_axis, mgf2_axis,
        )
        result = polarization_ray_trace(
            system, [0, 0, 1], [0, 0, 0], [1, 0],
            {"encodePropagationPhaseInP": True},
        )
        retardance[index] = jones_retardance_waves(coherent_final_jones(result))
    return retardance


def scalar_crossed_axis_retardance(
    wavelengths_um, quartz_thickness: float, mgf2_thickness: float
):
    wavelengths = np.asarray(wavelengths_um, dtype=float)
    n_oq, n_eq = quartz_optical_constants(wavelengths)
    n_om, n_em = mgf2_optical_constants(wavelengths)
    value = np.mod(
        ((n_em - n_om) * mgf2_thickness - (n_eq - n_oq) * quartz_thickness)
        / wavelengths,
        1,
    )
    return np.minimum(value, 1 - value)


def optimize_achromatic_wave_plate(
    wavelengths_um=None,
    target_waves: float = 0.25,
    initial_thicknesses=(230.537086, 190.217463),
    bounds=((1.0, 1.0), (2000.0, 2000.0)),
    gap_thickness: float = 444.5,
):
    """Optimize quartz and MgF2 thicknesses for broadband retardance."""
    if wavelengths_um is None:
        wavelengths_um = np.linspace(0.4, 0.7, 7)
    wavelengths = np.asarray(wavelengths_um, dtype=float)

    def residuals(thicknesses):
        return achromatic_retardance_spectrum(
            wavelengths, thicknesses[0], thicknesses[1],
            gap_thickness=gap_thickness,
        ) - target_waves

    return least_squares(
        residuals, initial_thicknesses, bounds=bounds,
        ftol=1e-12, xtol=1e-12, gtol=1e-12,
    )
