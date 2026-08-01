# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Sequential polarization ray tracing."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict

import numpy as np
from numpy.typing import ArrayLike

from .coordinates import double_pole_basis_vectors
from .geometry import intersect_surface
from .interfaces import (
    trace_isotropic_to_isotropic,
    trace_isotropic_to_uniaxial,
    trace_uniaxial_to_isotropic,
    trace_uniaxial_to_uniaxial,
)
from .linalg import make_k, normalize
from .models import (
    OpticalSystem,
    RayBranch,
    TraceOptions,
    TraceResult,
)


def _options(
    value: TraceOptions | Mapping[str, object] | None,
) -> TraceOptions:
    if value is None:
        return TraceOptions()
    if isinstance(value, TraceOptions):
        return value
    aliases = {
        "maxBranches": "max_branches",
        "minFlux": "min_flux",
        "minRelativeFlux": "min_relative_flux",
        "minAmplitude": "min_amplitude",
        "keepDiagnostics": "keep_diagnostics",
        "encodePropagationPhaseInP": "encode_propagation_phase_in_p",
        "maxInteractions": "max_interactions",
        "faceTolerance": "face_tolerance",
    }
    normalized = {aliases.get(key, key): item for key, item in value.items()}
    defaults = asdict(TraceOptions())
    defaults.update(normalized)
    return TraceOptions(**defaults)


def _initial_ray(
    system: OpticalSystem,
    k_incident: ArrayLike,
    position: ArrayLike,
    e_incident: ArrayLike,
) -> RayBranch:
    if not system.surfaces:
        raise ValueError("the optical system has no surfaces")
    k = normalize(k_incident)
    e_input = np.asarray(e_incident, dtype=np.complex128).reshape(-1)
    if e_input.size == 2:
        x_dp, y_dp = double_pole_basis_vectors(
            k, [0, 0, 1], [1, 0, 0]
        )
        field_e = e_input[0] * x_dp + e_input[1] * y_dp
        local_basis = {
            "x": x_dp,
            "y": y_dp,
            "basisDirection": k,
            "inputConvention": "doublePoleJones",
        }
    elif e_input.size == 3:
        field_e = e_input - np.dot(e_input, k) * k
        local_basis = {
            "x": None,
            "y": None,
            "basisDirection": k,
            "inputConvention": "global3D",
        }
    else:
        raise ValueError("E_in must contain two or three elements")
    if np.linalg.norm(field_e) == 0:
        raise ValueError("E_in must have nonzero transverse amplitude")

    first = system.surfaces[0]
    if first.material_type == "isotropic":
        n_incident = complex(first.index_data["n"])
        mode = "isotropic"
    else:
        n_incident = complex(first.index_data["nO"])
        mode = "input"
    mode_e = normalize(field_e)
    mode_h = n_incident * make_k(k) @ mode_e
    field_h = n_incident * make_k(k) @ field_e
    flux = float(np.real(np.dot(np.cross(field_e, np.conj(field_h)), k)))
    return RayBranch(
        id=1,
        mode=mode,
        branch_type="input",
        medium_type=first.material_type,
        position=np.asarray(position, dtype=np.complex128).reshape(3),
        k=k,
        s_direction=k,
        mode_e=mode_e,
        mode_h=mode_h,
        field_e=field_e,
        field_h=field_h,
        local_basis=local_basis,
        amplitude=float(np.linalg.norm(field_e)),
        flux=flux,
        active=True,
        metadata={"n": n_incident, "currentVertexZ": 0.0, "history": [1]},
    )


def _propagate(
    ray: RayBranch,
    hit: np.ndarray,
    wavelength: float,
    encode_phase: bool,
) -> None:
    segment = hit - ray.position
    phase_index = ray.metadata.get("phaseIndex", ray.metadata.get("n", 1.0))
    segment_opl = float(np.real(phase_index * np.dot(ray.k, segment)))
    ray.opl += segment_opl
    ray.position = np.asarray(hit, dtype=np.complex128)
    ray.metadata["lastSegmentOPL"] = segment_opl
    if not encode_phase:
        ray.metadata["lastSegmentPhase"] = None
        return
    phase = np.exp(1j * 2 * np.pi * segment_opl / wavelength)
    ray.field_e *= phase
    ray.field_h *= phase
    ray.metadata["lastSegmentPhase"] = phase
    required = {
        "P_interface",
        "P_beforePropagation",
        "propagationProjector",
    }
    if required.issubset(ray.metadata):
        projector = ray.metadata["propagationProjector"]
        p_propagated = (
            ray.metadata["P_interface"] - projector
        ) * phase + projector
        ray.p_matrix = p_propagated @ ray.metadata["P_beforePropagation"]
        ray.metadata["P_propagated"] = p_propagated


def _dispatch(
    ray: RayBranch,
    medium_in,
    medium_out,
    hit: np.ndarray,
    normal: np.ndarray,
    surface_index: int,
):
    case = (medium_in.material_type, medium_out.material_type)
    if case == ("isotropic", "isotropic"):
        return trace_isotropic_to_isotropic(
            ray, medium_in, medium_out, hit, normal, surface_index
        )
    if case == ("isotropic", "uniaxial"):
        return trace_isotropic_to_uniaxial(
            ray, medium_in, medium_out, hit, normal, surface_index
        )
    if case == ("uniaxial", "isotropic"):
        return trace_uniaxial_to_isotropic(
            ray, medium_in, medium_out, hit, normal, surface_index
        )
    if case == ("uniaxial", "uniaxial"):
        return trace_uniaxial_to_uniaxial(
            ray, medium_in, medium_out, hit, normal, surface_index
        )
    raise NotImplementedError(
        f"Python interface kernel {case[0]}To{case[1].title()} is pending"
    )


def polarization_ray_trace(
    system: OpticalSystem,
    k_incident: ArrayLike,
    position: ArrayLike,
    e_incident: ArrayLike,
    options: TraceOptions | Mapping[str, object] | None = None,
) -> TraceResult:
    """Trace a ray through currently supported planar interfaces."""

    trace_options = _options(options)
    result = TraceResult(system=system, options=trace_options)
    result.rays.append(
        _initial_ray(system, k_incident, position, e_incident)
    )
    active_ids = [1]
    for surface_index in range(len(system.surfaces) - 1):
        next_active_ids: list[int] = []
        for ray_id in active_ids:
            ray = result.rays[ray_id - 1]
            if not ray.active:
                continue
            surface = system.surfaces[surface_index]
            hit, normal, intercept_data = intersect_surface(surface, ray)
            _propagate(
                ray,
                hit,
                system.wavelength,
                trace_options.encode_propagation_phase_in_p,
            )
            interaction = _dispatch(
                ray,
                surface,
                system.surfaces[surface_index + 1],
                hit,
                normal,
                surface_index + 1,
            )
            interaction.diagnostics["interceptData"] = intercept_data
            result.interactions.append(interaction)
            parent_flux = max(abs(ray.flux), np.finfo(float).eps)
            target_z = float(ray.metadata.get("currentVertexZ", 0.0))
            target_z += surface.thickness
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
                child.surface_index = surface_index + 1
                child.metadata["currentVertexZ"] = target_z
                child.metadata["history"] = list(
                    ray.metadata.get("history", [ray_id])
                ) + [child.id]
                result.rays.append(child)
                is_last_interface = surface_index == len(system.surfaces) - 2
                if is_last_interface:
                    if child.branch_type == "transmitted":
                        next_active_ids.append(child.id)
                    else:
                        child.active = False
                    continue

                next_thickness = system.surfaces[surface_index + 1].thickness
                axial_direction = float(np.real(child.S[2]))
                routes_forward = (
                    abs(axial_direction) > 100 * np.finfo(float).eps
                    and next_thickness / axial_direction
                    > 100 * np.finfo(float).eps
                )
                if routes_forward:
                    next_active_ids.append(child.id)
                else:
                    child.active = False
            ray.active = False
        active_ids = next_active_ids
        if not active_ids:
            break
    result.final_ray_ids = active_ids
    return result
