# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Jones and retardance analysis helpers for traced ray trees."""

from __future__ import annotations

import numpy as np

from .coordinates import transform_p_to_jones
from .models import TraceResult


def retardance_by_qp_method(
    q_matrix: np.ndarray,
    p_matrix: np.ndarray,
    method: int = 0,
) -> float:
    """Calculate 3D retardance using Chipman's Q/P SVD construction."""

    q = np.asarray(q_matrix, dtype=np.complex128).reshape(3, 3)
    p = np.asarray(p_matrix, dtype=np.complex128).reshape(3, 3)
    u_matrix, _, vh_matrix = np.linalg.svd(np.linalg.solve(q, p))
    v_matrix = vh_matrix.conj().T
    eigenvalues = np.linalg.eigvals(np.linalg.solve(v_matrix, u_matrix))
    tolerance = 10 * np.finfo(float).eps
    candidates = np.flatnonzero(
        np.abs(np.abs(np.real(eigenvalues)) - 1) < tolerance
    )
    if candidates.size == 0:
        candidates = np.flatnonzero(
            np.abs(np.abs(np.real(eigenvalues)) - 1) < 1e-5
        )
    if candidates.size == 0:
        raise RuntimeError("could not identify the longitudinal eigenvalue")
    longitudinal = int(candidates[0])
    retarding = np.delete(eigenvalues, longitudinal)
    first, second = np.angle(retarding)
    high, low = (first, second) if first >= second else (second, first)
    if method == 0:
        return float(high - low)
    if method == 1:
        high_value = retarding[0] if first >= second else retarding[1]
        low_value = retarding[1] if first >= second else retarding[0]
        return float(np.angle(high_value - low_value))
    raise ValueError("method must be 0 or 1")


def jones_major_axis_ellipticity(field: np.ndarray):
    """Return major-axis orientation and unsigned ellipticity for Jones fields.

    ``field`` follows the MATLAB convention and must have shape ``(2, ...)``.
    Returned arrays have the trailing shape of ``field``; ``major_axis`` has
    shape ``(2, ...)``.
    """

    value = np.asarray(field, dtype=np.complex128)
    if value.ndim == 0 or value.shape[0] != 2:
        raise ValueError("field must have shape (2, ...)")
    e_x, e_y = value[0], value[1]
    a_x, a_y = np.abs(e_x), np.abs(e_y)
    cosine_phase = np.ones_like(a_x, dtype=float)
    nonzero = (a_x > 0) & (a_y > 0)
    cosine_phase[nonzero] = (
        np.real(e_x[nonzero] * np.conj(e_y[nonzero]))
        / (a_x[nonzero] * a_y[nonzero])
    )
    cosine_phase = np.clip(cosine_phase, -1.0, 1.0)
    psi = 0.5 * np.arctan2(
        2 * a_x * a_y * cosine_phase, a_x**2 - a_y**2
    )
    cosine_psi, sine_psi = np.cos(psi), np.sin(psi)
    cross_term = 2 * a_x * a_y * cosine_psi * sine_psi * cosine_phase
    major = np.sqrt(np.maximum(
        0.0, a_x**2 * cosine_psi**2 + a_y**2 * sine_psi**2 + cross_term
    ))
    minor = np.sqrt(np.maximum(
        0.0, a_x**2 * sine_psi**2 + a_y**2 * cosine_psi**2 - cross_term
    ))
    ellipticity = np.divide(
        minor, major, out=np.zeros_like(major), where=major > 0
    )
    major_axis = np.stack((cosine_psi, sine_psi), axis=0)
    return psi, ellipticity, major_axis, major, minor


def coherent_final_p(result: TraceResult) -> np.ndarray:
    """Coherently sum accumulated P matrices for all final branches."""
    return sum(
        (result.rays[ray_id - 1].P for ray_id in result.final_ray_ids),
        np.zeros((3, 3), dtype=np.complex128),
    )


def coherent_final_jones(
    result: TraceResult,
    k_in=(0, 0, 1),
    k_out=(0, 0, 1),
    coordinates=None,
    coordinates_out=None,
) -> np.ndarray:
    """Return a 2D Jones matrix from the coherent final P matrix."""
    if coordinates is None:
        coordinates = {
            "type": "doublePole", "a_loc": np.asarray(k_in),
            "x_o": np.array([1.0, 0.0, 0.0]),
        }
    return transform_p_to_jones(
        coherent_final_p(result), k_in, k_out, coordinates, coordinates_out
    )


def select_dominant_final_ray(result: TraceResult):
    """Return the final ray carrying the largest nonnegative flux."""

    if not result.final_ray_ids:
        raise ValueError("trace result contains no final rays")
    final_rays = [result.rays[ray_id - 1] for ray_id in result.final_ray_ids]
    ray = max(final_rays, key=lambda item: float(np.real(item.flux)))
    return ray, ray.id


def jones_retardance(jones: np.ndarray) -> float:
    """Return Chipman's principal retardance in radians.

    This implements PLAOS Equation 17.31 and remains valid when the Jones
    matrix contains diattenuation as well as retardance.
    """

    matrix = np.asarray(jones, dtype=np.complex128)
    determinant = np.linalg.det(matrix)
    if abs(determinant) < 1e-14:
        return np.nan
    numerator = abs(
        np.trace(matrix)
        + determinant / abs(determinant) * np.trace(matrix.conj().T)
    )
    denominator = 2 * np.sqrt(
        np.trace(matrix.conj().T @ matrix) + 2 * abs(determinant)
    )
    argument = float(np.clip(np.real(numerator / denominator), 0.0, 1.0))
    return float(2 * np.arccos(argument))


def jones_retardance_waves(jones: np.ndarray) -> float:
    """Return the principal retardance in the interval [0, 0.5] waves."""

    return jones_retardance(jones) / (2 * np.pi)
