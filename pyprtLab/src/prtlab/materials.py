# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Material eigenmode calculations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from .linalg import canonicalize_mode, make_k, normalize


@dataclass(slots=True)
class Mode:
    name: str
    q: np.ndarray
    k: np.ndarray
    e: np.ndarray
    h: np.ndarray
    s: np.ndarray
    phase_index: complex


def uniaxial_dielectric_tensor(
    n_o: float, n_e: float, optic_axis: ArrayLike
) -> np.ndarray:
    axis = normalize(optic_axis)
    identity = np.eye(3, dtype=np.complex128)
    return n_o**2 * identity + (n_e**2 - n_o**2) * np.outer(axis, axis)


def _mode_from_q(
    name: str,
    q: np.ndarray,
    epsilon: np.ndarray,
) -> Mode:
    matrix = epsilon + make_k(q) @ make_k(q)
    _, _, vh = np.linalg.svd(matrix)
    e = canonicalize_mode(vh.conj().T[:, -1])
    h = np.cross(q, e)
    poynting = np.real(np.cross(e, np.conj(h)))
    s = normalize(poynting)
    phase_index = np.sqrt(np.dot(q, q))
    k = q / phase_index
    return Mode(name, q, k, e, h, s, phase_index)


def _select_direction(
    candidates: list[Mode],
    normal: np.ndarray,
    direction_sign: int,
) -> Mode:
    valid = [
        mode
        for mode in candidates
        if direction_sign * np.real(np.dot(mode.s, normal)) > 0
    ]
    if not valid:
        raise RuntimeError("no uniaxial mode propagates in the requested direction")
    return max(
        valid,
        key=lambda mode: direction_sign * np.real(np.dot(mode.s, normal)),
    )


def uniaxial_modes_from_tangential_q(
    q_tangential: ArrayLike,
    normal: ArrayLike,
    n_o: float,
    n_e: float,
    optic_axis: ArrayLike,
    direction_sign: int = 1,
) -> tuple[Mode, Mode]:
    """Return ordinary and extraordinary modes sharing tangential q."""

    qt = np.asarray(q_tangential, dtype=np.complex128)
    ada = normalize(normal)
    axis = normalize(optic_axis)
    epsilon = uniaxial_dielectric_tensor(n_o, n_e, axis)

    ordinary_root = np.lib.scimath.sqrt(n_o**2 - np.dot(qt, qt))
    ordinary_candidates = [
        _mode_from_q("ordinary", qt + sign * ordinary_root * ada, epsilon)
        for sign in (1, -1)
    ]
    ordinary = _select_direction(ordinary_candidates, ada, direction_sign)

    axis_normal = np.dot(axis, ada)
    axis_tangent = np.dot(axis, qt)
    a = axis_normal**2 / n_o**2 + (1 - axis_normal**2) / n_e**2
    b = 2 * axis_tangent * axis_normal * (1 / n_o**2 - 1 / n_e**2)
    c = (
        axis_tangent**2 / n_o**2
        + (np.dot(qt, qt) - axis_tangent**2) / n_e**2
        - 1
    )
    extraordinary_roots = np.roots([a, b, c])
    extraordinary_candidates = [
        _mode_from_q("extraordinary", qt + root * ada, epsilon)
        for root in extraordinary_roots
    ]
    extraordinary = _select_direction(
        extraordinary_candidates, ada, direction_sign
    )
    return ordinary, extraordinary
