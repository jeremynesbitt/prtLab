# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Tests for the final public MATLAB helper ports."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from prtlab import (
    example_cell_phone_lens_ex9_system,
    plot_basis_vectors_on_sphere,
    plot_prt_lens_cross_section,
    plot_prt_ray_trace,
    polarization_ray_trace,
    retardance_by_qp_method,
    update_clear_apertures_from_ray_trace,
)


def _cell_phone_trace():
    system = example_cell_phone_lens_ex9_system()
    result = polarization_ray_trace(
        system,
        [0, 0, 1],
        [0, 0.741676704231623, -0.1],
        [0, 1],
    )
    return system, result


def test_qp_retardance_matches_matlab_diagonal_reference():
    q_matrix = np.eye(3)
    p_matrix = np.diag(np.exp(1j * np.array([0.2, 0.7, 0.0])))
    np.testing.assert_allclose(
        retardance_by_qp_method(q_matrix, p_matrix), 0.5, atol=1e-14
    )
    np.testing.assert_allclose(
        retardance_by_qp_method(q_matrix, p_matrix, method=1),
        2.0207963267948967,
        atol=1e-14,
    )


def test_clear_apertures_follow_intercepts_and_share_plane_pairs():
    system, result = _cell_phone_trace()
    _, apertures = update_clear_apertures_from_ray_trace(
        system, result, margin=1.1, minimum=0.2
    )
    expected = np.zeros(len(system.surfaces) - 1)
    for interaction in result.interactions:
        index = interaction.surface_index - 1
        expected[index] = max(
            expected[index], np.hypot(*np.real(interaction.position[:2]))
        )
    expected = np.maximum(1.1 * expected, 0.2)
    expected[8:10] = np.max(expected[8:10])
    np.testing.assert_allclose(apertures, expected)
    for index, aperture in enumerate(apertures):
        assert system.surfaces[index].surface_data["clearAperture"] == aperture


def test_2d_system_and_ray_plot_include_image_plane_and_final_ray():
    system, result = _cell_phone_trace()
    update_clear_apertures_from_ray_trace(system, result, margin=1.1)
    axis, system_handles = plot_prt_lens_cross_section(system)
    _, ray_handles = plot_prt_ray_trace(result, ax=axis, post_extend=0.5)
    final = result.rays[result.final_ray_ids[0] - 1]
    assert system_handles["elements"]
    assert ray_handles
    assert axis.get_xlim()[1] >= np.real(final.position[2])
    assert axis.get_xlabel() == "z (mm)"
    plt.close(axis.figure)


def test_basis_vector_sphere_plot_supports_one_or_two_bases():
    k_grid = np.array([[0, 0, 1], [0, 1, 0]], dtype=float)
    x_local = np.array([[1, 0, 0], [1, 0, 0]], dtype=float)
    y_local = np.cross(k_grid, x_local)
    axis, handles = plot_basis_vectors_on_sphere(k_grid, x_local, y_local)
    assert handles["sphere"] is not None
    assert handles["x"] is not None
    assert handles["y"] is not None
    np.testing.assert_allclose(axis.get_box_aspect(), axis.get_box_aspect()[0])
    plt.close(axis.figure)


def test_basis_vector_sphere_plot_rejects_empty_input():
    with pytest.raises(ValueError, match="at least one propagation vector"):
        plot_basis_vectors_on_sphere([], [], None)


def test_basis_vector_sphere_plot_accepts_matlab_view_direction():
    axis, _ = plot_basis_vectors_on_sphere(
        [[0, 0, 1]], [[1, 0, 0]], None,
        show_sphere=False, view=(0, 7, -90),
    )
    np.testing.assert_allclose(axis.azim, 90, atol=1e-12)
    np.testing.assert_allclose(
        axis.elev, np.rad2deg(np.arctan2(-90, 7)), atol=1e-12
    )
    plt.close(axis.figure)
