# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Sequential surface sag, intercept, and normal calculations."""

from __future__ import annotations

from typing import Any

import numpy as np

from .linalg import normalize
from .models import RayBranch, Surface


def odd_asphere_parameters(
    surface_data: dict[str, Any],
) -> tuple[float, np.ndarray]:
    if "kappa" in surface_data:
        kappa = float(surface_data["kappa"])
    elif "data" in surface_data:
        kappa = float(surface_data["data"][0])
    else:
        kappa = 0.0
    if "A" in surface_data:
        coefficients = surface_data["A"]
    elif "coefficients" in surface_data:
        coefficients = surface_data["coefficients"]
    elif "data" in surface_data:
        coefficients = surface_data["data"][1:]
    else:
        coefficients = []
    a = np.zeros(8, dtype=float)
    supplied = np.asarray(coefficients, dtype=float).reshape(-1)
    a[: min(8, supplied.size)] = supplied[:8]
    return kappa, a


def odd_asphere_sag(
    x: float,
    y: float,
    radius: float,
    kappa: float,
    coefficients: np.ndarray,
) -> float:
    radius_squared = x * x + y * y
    radial_coordinate = np.sqrt(radius_squared)
    if np.isinf(radius):
        base_sag = 0.0
    else:
        curvature = 1.0 / radius
        argument = max(1.0 - kappa * curvature**2 * radius_squared, 0.0)
        base_sag = (
            curvature * radius_squared / (1.0 + np.sqrt(argument))
        )
    polynomial_sag = sum(
        coefficients[power - 3] * radial_coordinate**power
        for power in range(3, 11)
    )
    return float(base_sag + polynomial_sag)


def surface_sag(x: float, y: float, surface: Surface) -> float:
    if surface.surface_type == "plane":
        return 0.0
    if surface.surface_type == "odd_asphere":
        kappa, coefficients = odd_asphere_parameters(surface.surface_data)
        return odd_asphere_sag(
            x, y, surface.radius, kappa, coefficients
        )
    raise NotImplementedError(
        f'surface type "{surface.surface_type}" is not implemented'
    )


def _sag_gradient(
    x: float,
    y: float,
    radius: float,
    kappa: float,
    coefficients: np.ndarray,
) -> tuple[float, float]:
    delta = 1e-7
    dx = (
        odd_asphere_sag(x + delta, y, radius, kappa, coefficients)
        - odd_asphere_sag(x - delta, y, radius, kappa, coefficients)
    ) / (2 * delta)
    dy = (
        odd_asphere_sag(x, y + delta, radius, kappa, coefficients)
        - odd_asphere_sag(x, y - delta, radius, kappa, coefficients)
    ) / (2 * delta)
    return dx, dy


def intersect_surface(
    surface: Surface,
    ray: RayBranch,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    previous_vertex = float(ray.metadata.get("currentVertexZ", 0.0))
    target_vertex = previous_vertex + surface.thickness
    direction = np.real_if_close(ray.S).astype(float)
    position = np.real_if_close(ray.position).astype(float)

    if surface.surface_type == "plane":
        if abs(direction[2]) <= 100 * np.finfo(float).eps:
            raise RuntimeError("ray is parallel to the next vertex plane")
        step = (target_vertex - position[2]) / direction[2]
        hit = position + step * direction
        normal = np.array([0.0, 0.0, 1.0])
        data = {
            "step": step,
            "targetZ": target_vertex,
            "surfaceType": surface.surface_type,
        }
    elif surface.surface_type == "odd_asphere":
        kappa, coefficients = odd_asphere_parameters(surface.surface_data)
        step = (
            (target_vertex - position[2]) / direction[2]
            if abs(direction[2]) > 1e-12
            else 1e-3
        )
        residual = np.nan
        for iterations in range(1, 301):
            point = position + step * direction
            sag = odd_asphere_sag(
                point[0], point[1], surface.radius, kappa, coefficients
            )
            residual = point[2] - target_vertex - sag
            d_sag_dx, d_sag_dy = _sag_gradient(
                point[0], point[1], surface.radius, kappa, coefficients
            )
            derivative = (
                direction[2]
                - d_sag_dx * direction[0]
                - d_sag_dy * direction[1]
            )
            if abs(derivative) < 1e-15:
                break
            step_delta = -residual / derivative
            step += step_delta
            if abs(step_delta) < 1e-12:
                break
        hit = position + step * direction
        d_sag_dx, d_sag_dy = _sag_gradient(
            hit[0], hit[1], surface.radius, kappa, coefficients
        )
        normal = normalize([-d_sag_dx, -d_sag_dy, 1.0]).real
        data = {
            "step": step,
            "targetZ": target_vertex,
            "surfaceType": surface.surface_type,
            "kappa": kappa,
            "A": coefficients,
            "iterations": iterations,
            "residual": residual,
        }
    else:
        raise NotImplementedError(
            f'surface type "{surface.surface_type}" is not implemented'
        )
    normal = normalize(normal).real
    if np.dot(normal, direction) < 0:
        normal = -normal
    return (
        np.asarray(hit, dtype=np.complex128),
        np.asarray(normal, dtype=np.complex128),
        data,
    )
