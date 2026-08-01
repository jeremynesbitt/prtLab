# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Isotropic-to-isotropic boundary-condition kernel."""

from __future__ import annotations

import numpy as np

from ..coordinates import calc_o, sp_basis
from ..linalg import make_k, normalize, snell_vector
from ..models import Interaction, RayBranch, Surface
from .common import boundary_matrix, build_child


def trace_isotropic_to_isotropic(
    incident: RayBranch,
    medium_in: Surface,
    medium_out: Surface,
    hit: np.ndarray,
    normal: np.ndarray,
    surface_index: int = 1,
) -> Interaction:
    ada = normalize(normal)
    k_in = normalize(incident.k)
    n1 = complex(incident.metadata.get("n", medium_in.index_data["n"]))
    n2 = complex(medium_out.index_data["n"])
    p_in, s_in = sp_basis(k_in, ada)
    h_in_s = n1 * make_k(k_in) @ s_in
    h_in_p = n1 * make_k(k_in) @ p_in

    k_out = snell_vector(k_in, ada, n1, n2)
    p_out, s_out = sp_basis(k_out, ada, s_in)
    h_out_s = n2 * make_k(k_out) @ s_out
    h_out_p = n2 * make_k(k_out) @ p_out

    k_ref = k_in - 2 * np.dot(k_in, ada) * ada
    p_ref, s_ref = sp_basis(k_ref, ada, s_in)
    h_ref_s = n1 * make_k(k_ref) @ s_ref
    h_ref_p = n1 * make_k(k_ref) @ p_ref

    s1 = s_in
    s2 = normalize(np.cross(ada, s1))
    f_matrix = boundary_matrix(
        [(s_out, h_out_s), (p_out, h_out_p)],
        [(s_ref, h_ref_s), (p_ref, h_ref_p)],
        s1,
        s2,
    )
    c_s = np.array(
        [np.dot(s1, s_in), np.dot(s2, s_in),
         np.dot(s1, h_in_s), np.dot(s2, h_in_s)]
    )
    c_p = np.array(
        [np.dot(s1, p_in), np.dot(s2, p_in),
         np.dot(s1, h_in_p), np.dot(s2, h_in_p)]
    )
    a_s = np.linalg.solve(f_matrix, c_s)
    a_p = np.linalg.solve(f_matrix, c_p)
    u_in = np.column_stack((s_in, p_in, incident.s_direction))
    p_transmitted = (
        np.column_stack(
            (
                a_s[0] * s_out + a_s[1] * p_out,
                a_p[0] * s_out + a_p[1] * p_out,
                k_out,
            )
        )
        @ np.linalg.inv(u_in)
    )
    s_ref_direction = normalize(np.real(np.cross(p_ref, np.conj(h_ref_p))))
    p_reflected = (
        np.column_stack(
            (
                a_s[2] * s_ref + a_s[3] * p_ref,
                a_p[2] * s_ref + a_p[3] * p_ref,
                s_ref_direction,
            )
        )
        @ np.linalg.inv(u_in)
    )
    o_in = calc_o(s_in, p_in, k_in)
    o_out = calc_o(s_out, p_out, k_out)
    o_ref = calc_o(s_ref, p_ref, k_ref)
    q_out = o_out @ np.linalg.inv(o_in)
    i_reflect = np.diag([1, -1, 1])
    q_ref = o_ref @ i_reflect @ np.linalg.inv(o_in)

    transmitted = build_child(
        incident, hit, medium_out, mode="isotropic",
        branch_type="transmitted", k=k_out, s_direction=k_out,
        p_matrix=p_transmitted, q_matrix=q_out, o_matrix=o_out,
        local_basis={"s": s_out, "p": p_out, "basisDirection": k_out},
        index=n2, flux_normal=ada,
        propagation_projector=np.outer(k_out, k_in),
        active=bool(
            np.max(np.abs(np.imag(k_out))) <= 100 * np.finfo(float).eps
            and np.max(np.abs(np.imag(n2))) <= 100 * np.finfo(float).eps
        ),
    )
    reflected = build_child(
        incident, hit, medium_in, mode="isotropic",
        branch_type="reflected", k=k_ref, s_direction=s_ref_direction,
        p_matrix=p_reflected, q_matrix=q_ref, o_matrix=o_ref,
        local_basis={"s": s_ref, "p": p_ref, "basisDirection": k_ref},
        index=n1, flux_normal=-ada,
        propagation_projector=np.outer(s_ref_direction, k_in),
    )
    return Interaction(
        surface_index=surface_index,
        case_name="isotropicToIsotropic",
        position=hit,
        normal=ada,
        incident=incident,
        children=[transmitted, reflected],
        frames={"Oin": o_in, "Oout": o_out},
        p_matrices={"transmitted": p_transmitted, "reflected": p_reflected},
        q_matrices={"transmitted": q_out, "reflected": q_ref},
        coefficients={"As": a_s, "Ap": a_p},
        diagnostics={
            "F": f_matrix,
            "Cs": c_s,
            "Cp": c_p,
            "boundaryResidual_s": f_matrix @ a_s - c_s,
            "boundaryResidual_p": f_matrix @ a_p - c_p,
        },
    )
