# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Isotropic-to-uniaxial boundary-condition kernel."""

from __future__ import annotations

import numpy as np

from ..coordinates import calc_o, sp_basis
from ..linalg import make_k, normalize
from ..materials import uniaxial_modes_from_tangential_q
from ..models import Interaction, RayBranch, Surface
from .common import boundary_matrix, build_child


def trace_isotropic_to_uniaxial(
    incident: RayBranch,
    medium_in: Surface,
    medium_out: Surface,
    hit: np.ndarray,
    normal: np.ndarray,
    surface_index: int = 1,
) -> Interaction:
    ada = normalize(normal)
    k_in = normalize(incident.k)
    n_inc = complex(incident.metadata.get("n", medium_in.index_data["n"]))
    n_o = float(medium_out.index_data["nO"])
    n_e = float(medium_out.index_data["nE"])
    optic_axis = normalize(medium_out.axis_data["opticAxis"])

    q_in = n_inc * k_in
    q_tangent = q_in - np.dot(q_in, ada) * ada
    ordinary, extraordinary = uniaxial_modes_from_tangential_q(
        q_tangent, ada, n_o, n_e, optic_axis, direction_sign=1
    )

    p_in, s_in = sp_basis(k_in, ada)
    h_in_s = n_inc * make_k(k_in) @ s_in
    h_in_p = n_inc * make_k(k_in) @ p_in
    k_ref = k_in - 2 * np.dot(k_in, ada) * ada
    p_ref, s_ref = sp_basis(k_ref, ada, s_in)
    h_ref_s = n_inc * make_k(k_ref) @ s_ref
    h_ref_p = n_inc * make_k(k_ref) @ p_ref
    s_ref_direction = normalize(np.real(np.cross(p_ref, np.conj(h_ref_p))))

    s1 = s_in
    s2 = normalize(np.cross(ada, s1))
    f_matrix = boundary_matrix(
        [(ordinary.e, ordinary.h), (extraordinary.e, extraordinary.h)],
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
    p_o = (
        np.column_stack(
            (a_s[0] * ordinary.e, a_p[0] * ordinary.e, ordinary.s)
        )
        @ np.linalg.inv(u_in)
    )
    p_e = (
        np.column_stack(
            (
                a_s[1] * extraordinary.e,
                a_p[1] * extraordinary.e,
                extraordinary.s,
            )
        )
        @ np.linalg.inv(u_in)
    )
    p_r = (
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
    p_o_basis, s_o_basis = sp_basis(ordinary.k, ada, s_in)
    o_out_o = calc_o(s_o_basis, p_o_basis, ordinary.k)
    p_e_basis, s_e_basis = sp_basis(extraordinary.s, ada, s_in)
    o_out_e = calc_o(s_e_basis, p_e_basis, extraordinary.s)
    p_r_basis, s_r_basis = sp_basis(k_ref, ada, s_in)
    o_ref = calc_o(s_r_basis, p_r_basis, k_ref)
    q_o = o_out_o @ np.linalg.inv(o_in)
    q_e = o_out_e @ np.linalg.inv(o_in)
    q_r = o_ref @ np.diag([1, -1, 1]) @ np.linalg.inv(o_in)

    child_o = build_child(
        incident, hit, medium_out, mode="ordinary",
        branch_type="transmitted", k=ordinary.k,
        s_direction=ordinary.s, p_matrix=p_o, q_matrix=q_o,
        o_matrix=o_out_o,
        local_basis={"s": s_o_basis, "p": p_o_basis,
                     "basisDirection": ordinary.k},
        index=ordinary.phase_index, mode_e=ordinary.e, mode_h=ordinary.h,
        flux_normal=ada, propagation_projector=np.outer(ordinary.k, k_in),
        metadata={"epsilon": None, "opticAxis": optic_axis},
    )
    child_e = build_child(
        incident, hit, medium_out, mode="extraordinary",
        branch_type="transmitted", k=extraordinary.k,
        s_direction=extraordinary.s, p_matrix=p_e, q_matrix=q_e,
        o_matrix=o_out_e,
        local_basis={"s": s_e_basis, "p": p_e_basis,
                     "basisDirection": extraordinary.s},
        index=extraordinary.phase_index, mode_e=extraordinary.e,
        mode_h=extraordinary.h, flux_normal=ada,
        phase_index=extraordinary.phase_index,
        propagation_projector=np.outer(extraordinary.s, k_in),
        metadata={
            "n_SE": np.dot(extraordinary.q, extraordinary.s),
            "epsilon": None,
            "opticAxis": optic_axis,
        },
    )
    child_r = build_child(
        incident, hit, medium_in, mode="isotropic",
        branch_type="reflected", k=k_ref, s_direction=s_ref_direction,
        p_matrix=p_r, q_matrix=q_r, o_matrix=o_ref,
        local_basis={"s": s_r_basis, "p": p_r_basis,
                     "basisDirection": k_ref},
        index=n_inc, flux_normal=-ada,
        propagation_projector=np.outer(s_ref_direction, k_in),
    )
    return Interaction(
        surface_index=surface_index,
        case_name="isotropicToUniaxial",
        position=hit,
        normal=ada,
        incident=incident,
        children=[child_o, child_e, child_r],
        frames={"Oin": o_in, "Oout_o": o_out_o, "Oout_e": o_out_e},
        p_matrices={"to_o": p_o, "to_e": p_e, "reflected": p_r},
        q_matrices={"ordinary": q_o, "extraordinary": q_e, "reflected": q_r},
        coefficients={
            "As": a_s, "Ap": a_p,
            "a_s_to": a_s[0], "a_s_te": a_s[1],
            "a_s_rs": a_s[2], "a_s_rp": a_s[3],
            "a_p_to": a_p[0], "a_p_te": a_p[1],
            "a_p_rs": a_p[2], "a_p_rp": a_p[3],
        },
        diagnostics={
            "F": f_matrix,
            "Cs": c_s,
            "Cp": c_p,
            "boundaryResidual_s": f_matrix @ a_s - c_s,
            "boundaryResidual_p": f_matrix @ a_p - c_p,
        },
    )


def trace_uniaxial_to_isotropic(
    incident: RayBranch,
    medium_in: Surface,
    medium_out: Surface,
    hit: np.ndarray,
    normal: np.ndarray,
    surface_index: int = 1,
) -> Interaction:
    """Trace one incident uniaxial eigenmode into an isotropic medium."""

    ada = normalize(normal)
    k_in = normalize(incident.k)
    s_in_direction = normalize(incident.s_direction)
    e_in = normalize(incident.mode_e)
    h_in = np.asarray(incident.mode_h, dtype=np.complex128)
    n_inc = complex(incident.metadata["n"])
    n_out = complex(medium_out.index_data["n"])
    n_o = float(medium_in.index_data["nO"])
    n_e = float(medium_in.index_data["nE"])
    optic_axis = normalize(medium_in.axis_data["opticAxis"])

    from ..linalg import snell_vector

    k_out = snell_vector(k_in, ada, n_inc, n_out)
    p_out, s_out = sp_basis(k_out, ada)
    h_out_s = n_out * make_k(k_out) @ s_out
    h_out_p = n_out * make_k(k_out) @ p_out

    q_in = n_inc * k_in
    q_tangent = q_in - np.dot(q_in, ada) * ada
    reflected_o, reflected_e = uniaxial_modes_from_tangential_q(
        q_tangent, ada, n_o, n_e, optic_axis, direction_sign=-1
    )

    interface_s = np.cross(k_in, ada)
    if np.linalg.norm(interface_s) <= 100 * np.finfo(float).eps:
        interface_s = np.array([1.0, 0.0, 0.0])
    interface_s = normalize(interface_s)
    interface_p = normalize(np.cross(ada, interface_s))
    f_matrix = boundary_matrix(
        [(s_out, h_out_s), (p_out, h_out_p)],
        [
            (reflected_o.e, reflected_o.h),
            (reflected_e.e, reflected_e.h),
        ],
        interface_s,
        interface_p,
    )
    c_mode = np.array(
        [
            np.dot(interface_s, e_in),
            np.dot(interface_p, e_in),
            np.dot(interface_s, h_in),
            np.dot(interface_p, h_in),
        ]
    )
    amplitudes = np.linalg.solve(f_matrix, c_mode)
    e_normal = normalize(np.cross(s_in_direction, e_in))
    modal_input = np.column_stack((e_in, e_normal, s_in_direction))

    e_transmitted = amplitudes[0] * s_out + amplitudes[1] * p_out
    h_transmitted = amplitudes[0] * h_out_s + amplitudes[1] * h_out_p
    p_transmitted = (
        np.column_stack((e_transmitted, np.zeros(3), k_out))
        @ np.linalg.inv(modal_input)
    )
    p_reflected_o = (
        np.column_stack(
            (
                amplitudes[2] * reflected_o.e,
                np.zeros(3),
                reflected_o.s,
            )
        )
        @ np.linalg.inv(modal_input)
    )
    p_reflected_e = (
        np.column_stack(
            (
                amplitudes[3] * reflected_e.e,
                np.zeros(3),
                reflected_e.s,
            )
        )
        @ np.linalg.inv(modal_input)
    )

    p_in_basis, s_in_basis = sp_basis(k_in, ada)
    o_in = calc_o(s_in_basis, p_in_basis, k_in)
    o_out = calc_o(s_out, p_out, k_out)
    p_ro_basis, s_ro_basis = sp_basis(reflected_o.k, ada)
    p_re_basis, s_re_basis = sp_basis(reflected_e.k, ada)
    o_ref_o = calc_o(s_ro_basis, p_ro_basis, reflected_o.k)
    o_ref_e = calc_o(s_re_basis, p_re_basis, reflected_e.k)
    q_out = o_out @ np.linalg.inv(o_in)
    reflection = np.diag([1, -1, 1])
    q_ref_o = o_ref_o @ reflection @ np.linalg.inv(o_in)
    q_ref_e = o_ref_e @ reflection @ np.linalg.inv(o_in)

    transmitted = build_child(
        incident, hit, medium_out, mode="isotropic",
        branch_type="transmitted", k=k_out, s_direction=k_out,
        p_matrix=p_transmitted, q_matrix=q_out, o_matrix=o_out,
        local_basis={"s": s_out, "p": p_out, "basisDirection": k_out},
        index=n_out, mode_e=e_transmitted, mode_h=h_transmitted,
        flux_normal=ada, propagation_projector=np.outer(k_out, k_in),
        metadata={"transmittedAmplitude": amplitudes},
    )
    child_o = build_child(
        incident, hit, medium_in, mode="ordinary",
        branch_type="reflected", k=reflected_o.k,
        s_direction=reflected_o.s, p_matrix=p_reflected_o,
        q_matrix=q_ref_o, o_matrix=o_ref_o,
        local_basis={"s": s_ro_basis, "p": p_ro_basis,
                     "basisDirection": reflected_o.k},
        index=reflected_o.phase_index, mode_e=reflected_o.e,
        mode_h=reflected_o.h, flux_normal=-ada,
        propagation_projector=np.outer(reflected_o.s, k_in),
        metadata={"reflectedAmplitude": amplitudes[2],
                  "opticAxis": optic_axis},
    )
    child_e = build_child(
        incident, hit, medium_in, mode="extraordinary",
        branch_type="reflected", k=reflected_e.k,
        s_direction=reflected_e.s, p_matrix=p_reflected_e,
        q_matrix=q_ref_e, o_matrix=o_ref_e,
        local_basis={"s": s_re_basis, "p": p_re_basis,
                     "basisDirection": reflected_e.k},
        index=reflected_e.phase_index, mode_e=reflected_e.e,
        mode_h=reflected_e.h, flux_normal=-ada,
        phase_index=reflected_e.phase_index,
        propagation_projector=np.outer(reflected_e.s, k_in),
        metadata={
            "reflectedAmplitude": amplitudes[3],
            "n_SE": np.dot(reflected_e.q, reflected_e.s),
            "opticAxis": optic_axis,
        },
    )
    return Interaction(
        surface_index=surface_index,
        case_name="uniaxialToIsotropic",
        position=hit,
        normal=ada,
        incident=incident,
        children=[transmitted, child_o, child_e],
        frames={"Oin": o_in, "Oout": o_out},
        p_matrices={
            "transmitted": p_transmitted,
            "reflectedOrdinary": p_reflected_o,
            "reflectedExtraordinary": p_reflected_e,
        },
        q_matrices={
            "transmitted": q_out,
            "reflectedOrdinary": q_ref_o,
            "reflectedExtraordinary": q_ref_e,
        },
        coefficients={
            "Am": amplitudes,
            "transmitted": amplitudes[:2],
            "reflectedOrdinary": amplitudes[2],
            "reflectedExtraordinary": amplitudes[3],
        },
        diagnostics={
            "F": f_matrix,
            "Cm": c_mode,
            "boundaryResidual_m": f_matrix @ amplitudes - c_mode,
        },
    )


def trace_uniaxial_to_uniaxial(
    incident: RayBranch,
    medium_in: Surface,
    medium_out: Surface,
    hit: np.ndarray,
    normal: np.ndarray,
    surface_index: int = 1,
) -> Interaction:
    """Trace one incident eigenmode between two uniaxial materials."""

    ada = normalize(normal)
    k_in = normalize(incident.k)
    s_in_direction = normalize(incident.s_direction)
    e_in = normalize(incident.mode_e)
    h_in = np.asarray(incident.mode_h, dtype=np.complex128)
    n_inc = complex(incident.metadata["n"])
    n_o_in = float(medium_in.index_data["nO"])
    n_e_in = float(medium_in.index_data["nE"])
    n_o_out = float(medium_out.index_data["nO"])
    n_e_out = float(medium_out.index_data["nE"])
    axis_in = normalize(medium_in.axis_data["opticAxis"])
    axis_out = normalize(medium_out.axis_data["opticAxis"])

    q_in = n_inc * k_in
    q_tangent = q_in - np.dot(q_in, ada) * ada
    transmitted_o, transmitted_e = uniaxial_modes_from_tangential_q(
        q_tangent, ada, n_o_out, n_e_out, axis_out, direction_sign=1
    )
    reflected_o, reflected_e = uniaxial_modes_from_tangential_q(
        q_tangent, ada, n_o_in, n_e_in, axis_in, direction_sign=-1
    )

    interface_s = np.cross(k_in, ada)
    if np.linalg.norm(interface_s) <= 100 * np.finfo(float).eps:
        interface_s = np.array([1.0, 0.0, 0.0])
    interface_s = normalize(interface_s)
    interface_p = normalize(np.cross(ada, interface_s))
    f_matrix = boundary_matrix(
        [
            (transmitted_o.e, transmitted_o.h),
            (transmitted_e.e, transmitted_e.h),
        ],
        [
            (reflected_o.e, reflected_o.h),
            (reflected_e.e, reflected_e.h),
        ],
        interface_s,
        interface_p,
    )
    c_mode = np.array(
        [
            np.dot(interface_s, e_in),
            np.dot(interface_p, e_in),
            np.dot(interface_s, h_in),
            np.dot(interface_p, h_in),
        ]
    )
    amplitudes = np.linalg.solve(f_matrix, c_mode)
    e_normal = normalize(np.cross(s_in_direction, e_in))
    modal_input = np.column_stack((e_in, e_normal, s_in_direction))

    def modal_p(amplitude, mode):
        return (
            np.column_stack(
                (amplitude * mode.e, np.zeros(3), mode.s)
            )
            @ np.linalg.inv(modal_input)
        )

    p_to = modal_p(amplitudes[0], transmitted_o)
    p_te = modal_p(amplitudes[1], transmitted_e)
    p_ro = modal_p(amplitudes[2], reflected_o)
    p_re = modal_p(amplitudes[3], reflected_e)

    p_in_basis, s_in_basis = sp_basis(k_in, ada)
    o_in = calc_o(s_in_basis, p_in_basis, k_in)

    def frame(mode, use_energy_direction=False):
        direction = mode.s if use_energy_direction else mode.k
        p_basis, s_basis = sp_basis(direction, ada, interface_s)
        return (
            calc_o(s_basis, p_basis, direction),
            p_basis,
            s_basis,
            direction,
        )

    o_to, p_to_basis, s_to_basis, direction_to = frame(transmitted_o)
    o_te, p_te_basis, s_te_basis, direction_te = frame(
        transmitted_e, use_energy_direction=True
    )
    o_ro, p_ro_basis, s_ro_basis, direction_ro = frame(reflected_o)
    o_re, p_re_basis, s_re_basis, direction_re = frame(reflected_e)
    q_to = o_to @ np.linalg.inv(o_in)
    q_te = o_te @ np.linalg.inv(o_in)
    q_ro = o_ro @ np.linalg.inv(o_in)
    q_re = o_re @ np.linalg.inv(o_in)

    child_to = build_child(
        incident, hit, medium_out, mode="ordinary",
        branch_type="transmitted", k=transmitted_o.k,
        s_direction=transmitted_o.s, p_matrix=p_to, q_matrix=q_to,
        o_matrix=o_to,
        local_basis={"s": s_to_basis, "p": p_to_basis,
                     "basisDirection": direction_to},
        index=transmitted_o.phase_index, mode_e=transmitted_o.e,
        mode_h=transmitted_o.h, flux_normal=ada,
        propagation_projector=np.outer(transmitted_o.k, k_in),
        metadata={"opticAxis": axis_out},
    )
    child_te = build_child(
        incident, hit, medium_out, mode="extraordinary",
        branch_type="transmitted", k=transmitted_e.k,
        s_direction=transmitted_e.s, p_matrix=p_te, q_matrix=q_te,
        o_matrix=o_te,
        local_basis={"s": s_te_basis, "p": p_te_basis,
                     "basisDirection": direction_te},
        index=transmitted_e.phase_index, mode_e=transmitted_e.e,
        mode_h=transmitted_e.h, flux_normal=ada,
        phase_index=transmitted_e.phase_index,
        propagation_projector=np.outer(transmitted_e.s, k_in),
        metadata={"n_SE": np.dot(transmitted_e.q, transmitted_e.s),
                  "opticAxis": axis_out},
    )
    child_ro = build_child(
        incident, hit, medium_in, mode="ordinary",
        branch_type="reflected", k=reflected_o.k,
        s_direction=reflected_o.s, p_matrix=p_ro, q_matrix=q_ro,
        o_matrix=o_ro,
        local_basis={"s": s_ro_basis, "p": p_ro_basis,
                     "basisDirection": direction_ro},
        index=reflected_o.phase_index, mode_e=reflected_o.e,
        mode_h=reflected_o.h, flux_normal=-ada,
        propagation_projector=np.outer(reflected_o.s, k_in),
        metadata={"opticAxis": axis_in},
    )
    child_re = build_child(
        incident, hit, medium_in, mode="extraordinary",
        branch_type="reflected", k=reflected_e.k,
        s_direction=reflected_e.s, p_matrix=p_re, q_matrix=q_re,
        o_matrix=o_re,
        local_basis={"s": s_re_basis, "p": p_re_basis,
                     "basisDirection": direction_re},
        index=reflected_e.phase_index, mode_e=reflected_e.e,
        mode_h=reflected_e.h, flux_normal=-ada,
        phase_index=reflected_e.phase_index,
        propagation_projector=np.outer(reflected_e.s, k_in),
        metadata={"n_SE": np.dot(reflected_e.q, reflected_e.s),
                  "opticAxis": axis_in},
    )
    return Interaction(
        surface_index=surface_index,
        case_name="uniaxialToUniaxial",
        position=hit,
        normal=ada,
        incident=incident,
        children=[child_to, child_te, child_ro, child_re],
        frames={"Oin": o_in, "Oout_o": o_to, "Oout_e": o_te},
        p_matrices={
            "to_o": p_to,
            "to_e": p_te,
            "reflectedOrdinary": p_ro,
            "reflectedExtraordinary": p_re,
        },
        q_matrices={
            "ordinary": q_to,
            "extraordinary": q_te,
            "reflectedOrdinary": q_ro,
            "reflectedExtraordinary": q_re,
        },
        coefficients={
            "Am": amplitudes,
            "transmittedOrdinary": amplitudes[0],
            "transmittedExtraordinary": amplitudes[1],
            "reflectedOrdinary": amplitudes[2],
            "reflectedExtraordinary": amplitudes[3],
        },
        diagnostics={
            "F": f_matrix,
            "Cm": c_mode,
            "boundaryResidual_m": f_matrix @ amplitudes - c_mode,
        },
    )
