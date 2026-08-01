# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Construction and transformation of extruded polygon optical solids."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from .linalg import normalize
from .models import Medium, Solid, SolidFace


def isotropic_medium(index: complex) -> Medium:
    return Medium("isotropic", {"n": index}, {})


def material_medium(
    index_data: complex | dict[str, Any], optic_axis: ArrayLike = (0, 1, 0)
) -> Medium:
    if np.isscalar(index_data):
        return isotropic_medium(complex(index_data))
    data = dict(index_data)
    if "n" in data:
        return Medium("isotropic", data, {})
    if "nO" in data and "nE" in data:
        return Medium(
            "uniaxial", data, {"opticAxis": normalize(optic_axis)}
        )
    raise ValueError("index_data must contain n, or both nO and nE")


def create_extruded_polygon_solid(
    name: str,
    vertices_yz: ArrayLike,
    face_names: list[str] | tuple[str, ...],
    width: float,
    outside: Medium,
    material: Medium,
    wavelength: float,
    wavelength_units: str = "um",
) -> Solid:
    vertices_yz = np.asarray(vertices_yz, dtype=float)
    if vertices_yz.ndim != 2 or vertices_yz.shape[0] != 2:
        raise ValueError("vertices_yz must have shape (2, N)")
    count = vertices_yz.shape[1]
    if count < 3 or len(face_names) != count:
        raise ValueError("provide at least three vertices and one name per edge")
    if width <= 0 or wavelength <= 0:
        raise ValueError("width and wavelength must be positive")

    center_yz = np.mean(vertices_yz, axis=1)
    center = np.array([0.0, *center_yz])
    faces: list[SolidFace] = []
    for index, face_name in enumerate(face_names):
        p0 = vertices_yz[:, index]
        p1 = vertices_yz[:, (index + 1) % count]
        vertices = np.array(
            [
                [-width / 2, p0[0], p0[1]],
                [width / 2, p0[0], p0[1]],
                [width / 2, p1[0], p1[1]],
                [-width / 2, p1[0], p1[1]],
            ],
            dtype=np.complex128,
        ).T
        normal = normalize(np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 1]))
        point = np.mean(vertices, axis=1)
        if np.real(np.dot(normal, point - center)) < 0:
            vertices = np.fliplr(vertices)
            normal = -normal
        faces.append(SolidFace(face_name, vertices, point, normal))
    return Solid(name, float(wavelength), wavelength_units, outside, material, faces)


def create_fresnel_rhomb_solid(
    index: float,
    *,
    length: float = 20,
    height: float = 8,
    width: float = 8,
    shear: float = 5,
    wavelength: float = 0.633,
    wavelength_units: str = "um",
) -> Solid:
    vertices = np.array([[0, height, height + shear, shear], [0, 0, length, length]])
    return create_extruded_polygon_solid(
        "Fresnel rhomb", vertices,
        ["entrance", "top", "exit", "bottom"], width,
        isotropic_medium(1.0), isotropic_medium(index), wavelength,
        wavelength_units,
    )


def create_right_triangular_prism_solid(
    index_data: complex | dict[str, Any],
    *,
    leg_y: float = 10,
    leg_z: float = 10,
    width: float = 8,
    origin: ArrayLike = (0, 0, 0),
    y_direction: int = 1,
    z_direction: int = 1,
    optic_axis: ArrayLike = (0, 1, 0),
    outside_index: float = 1.0,
    wavelength: float = 0.633,
    wavelength_units: str = "um",
    name: str = "Right triangular prism",
) -> Solid:
    origin = np.asarray(origin, dtype=float).reshape(3)
    a = origin[1:]
    b = a + [y_direction * leg_y, 0]
    c = a + [0, z_direction * leg_z]
    solid = create_extruded_polygon_solid(
        name, np.column_stack((a, b, c)),
        ["z leg", "hypotenuse", "y leg"], width,
        isotropic_medium(outside_index), material_medium(index_data, optic_axis),
        wavelength, wavelength_units,
    )
    solid.geometry = {
        "shape": "rightTriangle", "origin": origin, "legY": leg_y,
        "legZ": leg_z, "yDirection": y_direction,
        "zDirection": z_direction, "width": width,
    }
    return solid


def transform_solid(
    solid: Solid,
    rotation: ArrayLike | None = None,
    translation: ArrayLike = (0, 0, 0),
) -> Solid:
    transformed = deepcopy(solid)
    if rotation is None:
        rotation = np.eye(3)
    rotation = np.asarray(rotation, dtype=float).reshape(3, 3)
    translation = np.asarray(translation, dtype=float).reshape(3)
    if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-10) or not np.isclose(np.linalg.det(rotation), 1, atol=1e-10):
        raise ValueError("rotation must be a proper orthonormal matrix")
    for face in transformed.faces:
        face.vertices = rotation @ face.vertices + translation[:, None]
        face.point = rotation @ face.point + translation
        face.normal = rotation @ face.normal
    axis = transformed.material.axis_data.get("opticAxis")
    if axis is not None:
        transformed.material.axis_data["opticAxis"] = rotation @ axis
    if "origin" in transformed.geometry:
        transformed.geometry["origin"] = rotation @ transformed.geometry["origin"] + translation
    return transformed
