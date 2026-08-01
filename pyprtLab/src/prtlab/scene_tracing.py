# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Branched nearest-hit polarization tracing through multiple solids."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from numpy.typing import ArrayLike

from .models import OpticalSystem, RayBranch, Scene, TraceOptions, TraceResult
from .tracing import _dispatch, _initial_ray, _options, _propagate


def _point_in_face(point: np.ndarray, vertices: np.ndarray, normal: np.ndarray, tolerance: float) -> bool:
    for index in range(vertices.shape[1]):
        v1 = vertices[:, index]
        v2 = vertices[:, (index + 1) % vertices.shape[1]]
        side = np.dot(np.cross(v2 - v1, point - v1), normal)
        if np.real(side) < -tolerance:
            return False
    return True


def _solid_index(ray: RayBranch) -> int:
    return int(ray.metadata.get("solidIndex", 0))


def _intersect_scene(scene: Scene, ray: RayBranch, tolerance: float):
    current = _solid_index(ray)
    candidates = range(len(scene.solids)) if current == 0 else (current - 1,)
    best = None
    best_distance = np.inf
    for solid_index in candidates:
        solid = scene.solids[solid_index]
        for face_index, face in enumerate(solid.faces):
            denominator = np.dot(ray.k, face.normal)
            if current == 0:
                if np.real(denominator) >= -tolerance:
                    continue
            elif np.real(denominator) <= tolerance:
                continue
            distance = np.dot(face.point - ray.position, face.normal) / denominator
            if abs(np.imag(distance)) > tolerance or np.real(distance) <= tolerance or np.real(distance) >= best_distance:
                continue
            point = ray.position + distance * ray.k
            if _point_in_face(point, face.vertices, face.normal, tolerance):
                best_distance = float(np.real(distance))
                best = (point, solid_index, face_index, best_distance)
    return best


def _initial_scene_ray(scene: Scene, k_incident: ArrayLike, position: ArrayLike, e_incident: ArrayLike) -> RayBranch:
    holder = OpticalSystem(scene.wavelength)
    holder.surfaces.append(scene.outside.as_surface())
    ray = _initial_ray(holder, k_incident, position, e_incident)
    ray.metadata.update({"solidIndex": 0, "mediumName": "outside", "history": [1]})
    return ray


def polarization_ray_trace_scene(
    scene: Scene,
    k_incident: ArrayLike,
    position: ArrayLike,
    e_incident: ArrayLike,
    options: TraceOptions | Mapping[str, object] | None = None,
) -> TraceResult:
    """Trace all significant branches through disjoint homogeneous solids."""

    trace_options = _options(options)
    result = TraceResult(scene, trace_options)
    result.rays.append(_initial_scene_ray(scene, k_incident, position, e_incident))
    active_ids = [1]
    escaped_ids: list[int] = []

    for _ in range(trace_options.max_interactions):
        next_active: list[int] = []
        for ray_id in active_ids:
            ray = result.rays[ray_id - 1]
            if not ray.active:
                continue
            intersection = _intersect_scene(scene, ray, trace_options.face_tolerance)
            if intersection is None:
                ray.active = False
                escaped_ids.append(ray_id)
                continue
            hit, solid_index0, face_index0, distance = intersection
            solid = scene.solids[solid_index0]
            face = solid.faces[face_index0]
            current = _solid_index(ray)
            entering = current == 0
            medium_in = scene.outside if entering else solid.material
            medium_out = solid.material if entering else scene.outside
            normal = -face.normal if entering else face.normal
            output_solid_index = solid_index0 + 1 if entering else 0

            _propagate(ray, hit, scene.wavelength, trace_options.encode_propagation_phase_in_p)
            interaction = _dispatch(
                ray, medium_in.as_surface(face.coating_data),
                medium_out.as_surface(), hit, normal, face_index0 + 1,
            )
            interaction.diagnostics["interceptData"] = {
                "solidIndex": solid_index0 + 1,
                "solidName": solid.name,
                "faceIndex": face_index0 + 1,
                "faceName": face.name,
                "distance": distance,
                "direction": "entering" if entering else "exiting",
            }
            result.interactions.append(interaction)
            parent_flux = max(abs(ray.flux), np.finfo(float).eps)
            parent_history = list(ray.metadata.get("history", [ray_id]))
            for child in interaction.children:
                relative_flux = abs(child.flux) / parent_flux
                if (
                    not child.active
                    or abs(child.flux) <= trace_options.min_flux
                    or child.amplitude <= trace_options.min_amplitude
                    or relative_flux < trace_options.min_relative_flux
                    or len(result.rays) >= trace_options.max_branches
                ):
                    continue
                child.id = len(result.rays) + 1
                child.parent_id = ray_id
                child.surface_index = face_index0 + 1
                child_solid_index = current if child.branch_type == "reflected" else output_solid_index
                child.metadata.update(
                    {
                        "solidIndex": child_solid_index,
                        "mediumName": "outside" if child_solid_index == 0 else solid.name,
                        "solidFaceIndex": face_index0 + 1,
                        "solidFaceName": face.name,
                        "sceneSolidIndex": solid_index0 + 1,
                        "sceneSolidName": solid.name,
                        "history": parent_history + [child.id],
                    }
                )
                result.rays.append(child)
                next_active.append(child.id)
            ray.active = False
        active_ids = next_active
        if not active_ids:
            break
    result.final_ray_ids = list(dict.fromkeys(escaped_ids + active_ids))
    return result


def polarization_ray_trace_solid(
    solid,
    k_incident: ArrayLike,
    position: ArrayLike,
    e_incident: ArrayLike,
    options: TraceOptions | Mapping[str, object] | None = None,
) -> TraceResult:
    """Compatibility wrapper for tracing a scene containing one solid."""

    from .scenes import add_solid, create_scene

    scene = create_scene(
        solid.wavelength, solid.outside,
        wavelength_units=solid.wavelength_units, name=solid.name,
    )
    add_solid(scene, solid)
    return polarization_ray_trace_scene(scene, k_incident, position, e_incident, options)
