# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Polarization coordinate-system utilities."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from .linalg import as_vector3, normalize


def sp_basis(
    direction: ArrayLike,
    normal: ArrayLike,
    s_hint: ArrayLike | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    k = normalize(direction)
    ada = normalize(normal)
    s = np.cross(k, ada)
    if np.linalg.norm(s) <= 100 * np.finfo(float).eps:
        if s_hint is not None:
            s = as_vector3(s_hint)
            s = s - np.dot(s, k) * k
        else:
            trial = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(trial, k)) > 0.9:
                trial = np.array([0.0, 1.0, 0.0])
            s = trial - np.dot(trial, k) * k
    s = normalize(s)
    p = normalize(np.cross(k, s))
    return p, s


def calc_o(s: ArrayLike, p: ArrayLike, k: ArrayLike) -> np.ndarray:
    return np.column_stack((normalize(s), normalize(p), normalize(k)))


def double_pole_basis_vectors(
    k: ArrayLike,
    a_loc: ArrayLike,
    x_o: ArrayLike,
) -> tuple[np.ndarray, np.ndarray]:
    direction = normalize(k)
    pole = normalize(a_loc)
    x_reference = normalize(x_o)
    if abs(np.dot(pole, x_reference)) >= 1e-6:
        raise ValueError("x_o must be perpendicular to a_loc")
    y_reference = np.cross(pole, x_reference)
    cosine = float(np.real(np.dot(direction, pole)))
    if cosine > 1 - 1e-10 or cosine < -1 + 1e-10:
        return x_reference, y_reference
    rotation_axis = normalize(np.cross(direction, pole))
    theta = -np.arccos(np.clip(cosine, -1, 1))

    def rodrigues(vector: np.ndarray) -> np.ndarray:
        return (
            vector * np.cos(theta)
            + np.cross(rotation_axis, vector) * np.sin(theta)
            + rotation_axis
            * np.dot(rotation_axis, vector)
            * (1 - np.cos(theta))
        )

    return normalize(rodrigues(x_reference)), normalize(
        rodrigues(y_reference)
    )


def dipole_basis_vectors(
    k: ArrayLike,
    a_loc: ArrayLike,
) -> tuple[np.ndarray, np.ndarray]:
    """Return Chipman's dipole basis vectors for one propagation direction."""

    direction = normalize(k)
    pole = normalize(a_loc)
    axis_cross_k = np.cross(pole, direction)
    magnitude = np.linalg.norm(axis_cross_k)
    if magnitude > 1e-6:
        x_local = axis_cross_k / magnitude
        return x_local, normalize(np.cross(direction, x_local))

    reference = (
        np.array([1.0, 0.0, 0.0])
        if abs(direction[0]) < 0.9
        else np.array([0.0, 1.0, 0.0])
    )
    x_local = normalize(np.cross(direction, reference))
    return x_local, normalize(np.cross(direction, x_local))


def _jones_local_basis(
    direction: ArrayLike,
    coordinates: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    if "type" not in coordinates:
        raise ValueError("coordinate mapping must contain a 'type' field")
    coordinate_type = str(coordinates["type"]).lower()
    if coordinate_type == "doublepole":
        _require_coordinate_fields(coordinates, coordinate_type, "a_loc", "x_o")
        return double_pole_basis_vectors(
            direction, coordinates["a_loc"], coordinates["x_o"]
        )
    if coordinate_type == "dipole":
        _require_coordinate_fields(coordinates, coordinate_type, "a_loc")
        return dipole_basis_vectors(direction, coordinates["a_loc"])
    if coordinate_type == "sp":
        _require_coordinate_fields(coordinates, coordinate_type, "normal")
        k = normalize(direction)
        normal = normalize(coordinates["normal"])
        s_local = np.cross(k, normal)
        if np.linalg.norm(s_local) < 1e-12:
            raise ValueError(
                "s/p coordinates are singular because k is parallel "
                "to the surface normal"
            )
        s_local = normalize(s_local)
        return s_local, normalize(np.cross(k, s_local))
    raise ValueError(
        f"unknown coordinate type {coordinates['type']!r}; "
        "use 'doublePole', 'dipole', or 'sp'"
    )


def _require_coordinate_fields(
    coordinates: dict[str, Any], coordinate_type: str, *field_names: str
) -> None:
    for field_name in field_names:
        if field_name not in coordinates:
            raise ValueError(
                f"coordinate type {coordinate_type!r} requires "
                f"field {field_name!r}"
            )
        value = np.asarray(coordinates[field_name]).reshape(-1)
        if value.size != 3 or np.linalg.norm(value) == 0:
            raise ValueError(
                f"coordinate field {field_name!r} must be a nonzero 3-vector"
            )


def transform_p_to_jones(
    p_matrix: ArrayLike,
    k_in: ArrayLike,
    k_out: ArrayLike,
    coordinates: dict[str, Any],
    coordinates_out: dict[str, Any] | None = None,
) -> np.ndarray:
    """Convert a 3D polarization ray-tracing matrix to a Jones matrix.

    ``coordinates`` and ``coordinates_out`` accept ``doublePole``, ``dipole``,
    or ``sp`` coordinate mappings. When the output mapping is omitted, the
    input mapping is used at both pupils.
    """

    output_coordinates = coordinates if coordinates_out is None else coordinates_out
    normalized_k_in = normalize(k_in)
    normalized_k_out = normalize(k_out)
    x_in, y_in = _jones_local_basis(normalized_k_in, coordinates)
    x_out, y_out = _jones_local_basis(normalized_k_out, output_coordinates)
    u_in = np.column_stack((x_in, y_in, normalized_k_in))
    u_out = np.column_stack((x_out, y_out, normalized_k_out))
    transformed = np.linalg.solve(
        u_out, np.asarray(p_matrix, dtype=np.complex128) @ u_in
    )
    return transformed[:2, :2]
