# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Shared interface and child-ray construction."""

from __future__ import annotations

from typing import Any

import numpy as np

from ..linalg import make_k, normalize
from ..models import RayBranch, Surface


def material_index(surface: Surface, mode: str = "") -> complex:
    if surface.material_type == "isotropic":
        return complex(surface.index_data["n"])
    if mode == "extraordinary":
        return complex(surface.index_data["nE"])
    return complex(surface.index_data["nO"])


def build_child(
    parent: RayBranch,
    hit: np.ndarray,
    medium: Surface,
    *,
    mode: str,
    branch_type: str,
    k: np.ndarray,
    s_direction: np.ndarray,
    p_matrix: np.ndarray,
    q_matrix: np.ndarray,
    o_matrix: np.ndarray,
    local_basis: dict[str, Any],
    index: complex,
    mode_e: np.ndarray | None = None,
    mode_h: np.ndarray | None = None,
    flux_normal: np.ndarray,
    phase_index: complex | None = None,
    propagation_projector: np.ndarray | None = None,
    metadata: dict[str, Any] | None = None,
    active: bool = True,
) -> RayBranch:
    field_e = p_matrix @ parent.field_e
    if mode_e is None:
        resolved_e = (
            normalize(field_e)
            if np.linalg.norm(field_e) > 100 * np.finfo(float).eps
            else np.zeros(3)
        )
    else:
        resolved_e = np.asarray(mode_e, dtype=np.complex128)
    if mode_h is None:
        resolved_h = index * make_k(k) @ resolved_e
    else:
        resolved_h = np.asarray(mode_h, dtype=np.complex128)
    denominator = np.vdot(resolved_e, resolved_e)
    modal_scale = (
        0.0 if abs(denominator) < np.finfo(float).eps
        else np.vdot(resolved_e, field_e) / denominator
    )
    field_h = modal_scale * resolved_h
    flux = float(
        np.real(np.dot(np.cross(field_e, np.conj(field_h)), flux_normal))
    )
    child_metadata = {} if metadata is None else dict(metadata)
    child_metadata.update(
        {
            "n": index,
            "phaseIndex": index if phase_index is None else phase_index,
            "modalScale": modal_scale,
            "P_interface": p_matrix,
            "P_beforePropagation": parent.p_matrix,
            "IoutField": flux,
        }
    )
    if propagation_projector is not None:
        child_metadata["propagationProjector"] = propagation_projector
    return RayBranch(
        mode=mode,
        branch_type=branch_type,
        medium_type=medium.material_type,
        position=np.asarray(hit, dtype=np.complex128),
        k=np.asarray(k, dtype=np.complex128),
        s_direction=np.asarray(s_direction, dtype=np.complex128),
        mode_e=resolved_e,
        mode_h=resolved_h,
        field_e=field_e,
        field_h=field_h,
        p_matrix=p_matrix @ parent.p_matrix,
        q_matrix=q_matrix @ parent.q_matrix,
        o_matrix=o_matrix,
        local_basis=local_basis,
        amplitude=float(np.linalg.norm(field_e)),
        flux=flux,
        opl=parent.opl,
        active=active,
        metadata=child_metadata,
    )


def boundary_matrix(
    transmitted: list[tuple[np.ndarray, np.ndarray]],
    reflected: list[tuple[np.ndarray, np.ndarray]],
    s1: np.ndarray,
    s2: np.ndarray,
) -> np.ndarray:
    columns = []
    for e, h in transmitted:
        columns.append(
            np.array([np.dot(s1, e), np.dot(s2, e),
                      np.dot(s1, h), np.dot(s2, h)])
        )
    for e, h in reflected:
        columns.append(
            -np.array([np.dot(s1, e), np.dot(s2, e),
                       np.dot(s1, h), np.dot(s2, h)])
        )
    return np.column_stack(columns)
