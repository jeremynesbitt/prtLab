# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Optical-system construction helpers."""

from __future__ import annotations

from typing import Any

import numpy as np

from .models import OpticalSystem, Surface


def create_optical_system(wavelength: float) -> OpticalSystem:
    if wavelength <= 0:
        raise ValueError("wavelength must be positive")
    return OpticalSystem(wavelength=float(wavelength))


def add_surface(
    system: OpticalSystem,
    radius: float,
    thickness: float,
    surface_type: str,
    surface_data: dict[str, Any],
    material_type: str,
    index_data: dict[str, Any],
    axis_data: dict[str, Any],
    coating_data: dict[str, Any] | None = None,
) -> OpticalSystem:
    system.surfaces.append(
        Surface(
            radius=float(radius),
            thickness=float(thickness),
            surface_type=str(surface_type),
            surface_data=dict(surface_data),
            material_type=str(material_type),
            index_data=dict(index_data),
            axis_data=dict(axis_data),
            coating_data=(
                {"type": "none"}
                if coating_data is None
                else dict(coating_data)
            ),
        )
    )
    return system


def update_clear_apertures_from_ray_trace(
    system: OpticalSystem,
    ray_outputs,
    *,
    margin: float = 1.0,
    minimum: float = 0.0,
    surfaces=None,
) -> tuple[OpticalSystem, np.ndarray]:
    """Set plotting apertures from maximum radial ray intercepts."""

    if margin < 0 or minimum < 0:
        raise ValueError("margin and minimum must be nonnegative")
    outputs = (
        list(ray_outputs)
        if isinstance(ray_outputs, (list, tuple))
        else [ray_outputs]
    )
    count = max(0, len(system.surfaces) - 1)
    apertures = np.zeros(count, dtype=float)
    for result in outputs:
        for interaction in result.interactions:
            index = interaction.surface_index - 1
            if 0 <= index < count:
                hit = np.real(interaction.position)
                apertures[index] = max(apertures[index], np.hypot(hit[0], hit[1]))
    apertures = np.maximum(margin * apertures, minimum)

    index = 0
    while index < count - 1:
        if (
            system.surfaces[index].surface_type == "plane"
            and system.surfaces[index + 1].surface_type == "plane"
        ):
            pair = max(apertures[index], apertures[index + 1])
            apertures[index:index + 2] = pair
            index += 2
        else:
            index += 1

    selected = range(count) if surfaces is None else surfaces
    for surface_index in selected:
        if not 0 <= surface_index < count:
            continue
        system.surfaces[surface_index].surface_data["clearAperture"] = float(
            apertures[surface_index]
        )
    return system, apertures
