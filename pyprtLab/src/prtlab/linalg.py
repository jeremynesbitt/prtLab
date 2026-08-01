# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Linear-algebra helpers with optical complex-vector conventions."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

Array = NDArray[np.complex128]


def as_vector3(value: ArrayLike) -> Array:
    vector = np.asarray(value, dtype=np.complex128).reshape(3)
    return vector


def prt_norm(vector: ArrayLike) -> complex:
    """Algebraic vector length, without complex conjugation."""

    value = np.asarray(vector, dtype=np.complex128)
    return np.sqrt(np.sum(value * value))


def normalize(vector: ArrayLike, *, algebraic: bool = False) -> Array:
    value = np.asarray(vector, dtype=np.complex128)
    length = prt_norm(value) if algebraic else np.linalg.norm(value)
    if abs(length) <= 100 * np.finfo(float).eps:
        raise ValueError("cannot normalize a zero-length vector")
    return value / length


def make_k(vector: ArrayLike) -> Array:
    x, y, z = as_vector3(vector)
    return np.array(
        [[0, -z, y], [z, 0, -x], [-y, x, 0]],
        dtype=np.complex128,
    )


def calc_k_from_theta_x_theta_y(
    theta_x_deg: float, theta_y_deg: float
) -> Array:
    theta_x = np.deg2rad(theta_x_deg)
    theta_y = np.deg2rad(theta_y_deg)
    vector = np.array(
        [np.tan(theta_x), np.tan(theta_y), 1.0],
        dtype=np.complex128,
    )
    return normalize(vector)


def get_theta_x_theta_y_from_k(vector: ArrayLike) -> tuple[float, float]:
    """Recover Chipman's x/y direction angles in radians."""

    k = normalize(vector)
    if np.isclose(k[2], 0.0):
        raise ValueError("theta_x and theta_y are undefined when k[2] is zero")
    theta_x = np.arctan(np.real(k[0] / k[2]))
    theta_y = np.arctan(np.real(k[1] / k[2]))
    return float(theta_x), float(theta_y)


def snell_vector(
    k_incident: ArrayLike,
    normal: ArrayLike,
    n_incident: complex,
    n_out: complex,
) -> Array:
    """Vector Snell law, including complex evanescent wavevectors."""

    k_in = normalize(k_incident, algebraic=np.iscomplexobj(k_incident))
    ada = normalize(normal)
    eta = n_incident / n_out
    tangential = eta * (k_in - np.dot(k_in, ada) * ada)
    normal_component = np.lib.scimath.sqrt(1 - np.dot(tangential, tangential))
    if np.real(np.dot(k_in, ada)) < 0:
        normal_component = -normal_component
    return tangential + normal_component * ada


def canonicalize_mode(vector: ArrayLike) -> Array:
    """Choose a deterministic sign/phase for a modal eigenvector."""

    value = normalize(vector)
    pivot = int(np.argmax(np.abs(value)))
    if abs(value[pivot]) > 0:
        value = value * np.exp(-1j * np.angle(value[pivot]))
    if np.real(value[pivot]) < 0:
        value = -value
    return np.real_if_close(value).astype(np.complex128)
