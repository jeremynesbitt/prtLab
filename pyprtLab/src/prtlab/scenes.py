# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Multiple-solid optical scene construction."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from .models import Medium, Scene, Solid
from .solids import create_right_triangular_prism_solid, isotropic_medium


def create_scene(
    wavelength: float,
    outside_medium: float | Medium = 1.0,
    *,
    wavelength_units: str = "um",
    name: str = "prtLab scene",
) -> Scene:
    outside = isotropic_medium(outside_medium) if np.isscalar(outside_medium) else outside_medium
    return Scene(float(wavelength), outside, wavelength_units, name)


def add_solid(scene: Scene, solid: Solid) -> Scene:
    if not np.isclose(solid.wavelength, scene.wavelength, atol=100 * np.finfo(float).eps * max(solid.wavelength, scene.wavelength)):
        raise ValueError("all scene solids must use the scene wavelength")
    if solid.wavelength_units != scene.wavelength_units:
        raise ValueError("all scene solids must use the scene length units")
    solid.outside = scene.outside
    scene.solids.append(solid)
    return scene


def create_glan_taylor_polarizer_scene(
    index_data: dict[str, Any],
    *,
    leg_y: float = 10,
    leg_z: float = 10,
    width: float = 8,
    air_gap: float = 0.02,
    origin: ArrayLike = (0, 0, 0),
    optic_axis: ArrayLike = (1, 0, 0),
    outside_index: float = 1.0,
    wavelength: float = 0.633,
    wavelength_units: str = "um",
    name: str = "Glan-Taylor polarizer",
) -> Scene:
    origin = np.asarray(origin, dtype=float).reshape(3)
    scene = create_scene(wavelength, outside_index, wavelength_units=wavelength_units, name=name)
    prism1 = create_right_triangular_prism_solid(
        index_data, leg_y=leg_y, leg_z=leg_z, width=width, origin=origin,
        optic_axis=optic_axis, outside_index=outside_index,
        wavelength=wavelength, wavelength_units=wavelength_units,
        name="Glan-Taylor prism 1",
    )
    hypotenuse = next(face for face in prism1.faces if face.name == "hypotenuse")
    prism2_origin = origin + [0, leg_y, leg_z] + air_gap * np.real(hypotenuse.normal)
    prism2 = create_right_triangular_prism_solid(
        index_data, leg_y=leg_y, leg_z=leg_z, width=width,
        origin=prism2_origin, y_direction=-1, z_direction=-1,
        optic_axis=optic_axis, outside_index=outside_index,
        wavelength=wavelength, wavelength_units=wavelength_units,
        name="Glan-Taylor prism 2",
    )
    add_solid(scene, prism1)
    add_solid(scene, prism2)
    scene.geometry = {
        "type": "GlanTaylor", "airGap": air_gap, "legY": leg_y,
        "legZ": leg_z, "width": width,
        "opticAxis": np.asarray(optic_axis, dtype=float) / np.linalg.norm(optic_axis),
    }
    return scene
