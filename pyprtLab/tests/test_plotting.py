# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from prtlab import (
    create_fresnel_rhomb_solid,
    plot_jones_pupil,
    plot_jones_vector,
    plot_polarization_ellipses_across_field,
    plot_polarization_glyph_3d,
    plot_prt_system_3d,
    polarization_ellipse_points_2d,
    polarization_ray_trace_solid,
    prepare_jones_pupil,
)


def _rhomb_result():
    solid = create_fresnel_rhomb_solid(
        1.495, length=14, height=8, width=8, shear=11.02085,
        wavelength=0.589,
    )
    result = polarization_ray_trace_solid(
        solid,
        [0, np.sin(np.deg2rad(0.1)), np.cos(np.deg2rad(0.1))],
        [0, 5, -2], [1, 1] / np.sqrt(2),
        {"maxInteractions": 8, "minAmplitude": 1e-4},
    )
    return solid, result


def test_circular_polarization_glyph_has_equal_axes():
    figure = plt.figure()
    axis = figure.add_subplot(projection="3d")
    artists = plot_polarization_glyph_3d(
        axis, [0, 0, 0], [1, 1j, 0] / np.sqrt(2), samples=201
    )
    x, y, z = artists[0].get_data_3d()
    np.testing.assert_allclose(np.ptp(x), np.ptp(y), rtol=2e-4)
    np.testing.assert_allclose(z, 0, atol=1e-14)
    plt.close(figure)


def test_plot_prt_system_3d_creates_geometry_rays_and_glyphs():
    solid, result = _rhomb_result()
    axis, handles = plot_prt_system_3d(
        solid, result,
        {"RaySelection": "dominant", "PostExtend": 2, "ShowLegend": True},
    )
    assert handles["surfaces"]
    assert handles["rays"]
    assert handles["glyphs"]
    assert axis.get_xlabel() == "x (um)"
    plt.close(axis.figure)


def test_sequential_trace_records_plot_history():
    from prtlab import add_surface, create_optical_system, polarization_ray_trace

    system = create_optical_system(0.55)
    add_surface(system, np.inf, 1, "plane", {}, "isotropic", {"n": 1.0}, {})
    add_surface(system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.5}, {})
    result = polarization_ray_trace(system, [0, 0, 1], [0, 0, 0], [1, 0])
    for ray_id in result.final_ray_ids:
        history = result.rays[ray_id - 1].metadata["history"]
        assert history[0] == 1
        assert history[-1] == ray_id


def test_2d_circular_polarization_ellipse_has_equal_axes():
    points = polarization_ellipse_points_2d(
        np.array([1, 1j]) / np.sqrt(2), samples=201
    )
    np.testing.assert_allclose(np.ptp(points[0]), np.ptp(points[1]), rtol=2e-4)


def test_plot_jones_vector_creates_ellipse_and_arrow():
    axis, artists = plot_jones_vector([1, 0.4j])
    assert len(artists) == 2
    assert axis.get_aspect() == 1.0
    plt.close(axis.figure)


def test_prepare_jones_pupil_absorbs_real_sign_into_amplitude():
    jones = np.zeros((3, 3, 2, 2), dtype=complex)
    pattern = np.array([[-1, -1, 1], [-1, 1, 1], [-1, 1, 1]], dtype=float)
    jones[:, :, 0, 1] = pattern * np.exp(0.2j)
    jones[:, :, 0, 0] = 1
    jones[:, :, 1, 1] = 1
    mask = np.ones((3, 3), dtype=bool)
    data = prepare_jones_pupil(jones, mask)
    np.testing.assert_allclose(
        np.abs(data["signed_amplitude"][:, :, 0, 1]), 1, atol=1e-14
    )
    valid_phase = data["phase"][:, :, 0, 1]
    np.testing.assert_allclose(valid_phase, 0, atol=1e-14)
    assert set(np.unique(np.sign(data["signed_amplitude"][:, :, 0, 1]))) == {-1.0, 1.0}


def test_plot_jones_pupil_returns_two_matrix_figures():
    jones = np.zeros((5, 5, 2, 2), dtype=complex)
    jones[:, :, 0, 0] = 1
    jones[:, :, 1, 1] = np.exp(0.3j)
    amplitude, phase, data = plot_jones_pupil(jones)
    assert amplitude[1].shape == (2, 2)
    assert phase[1].shape == (2, 2)
    assert data["pupil_mask"].shape == (5, 5)
    plt.close(amplitude[0])
    plt.close(phase[0])


def test_plot_polarization_ellipses_across_field_counts_masked_samples():
    coordinates = np.linspace(-1, 1, 3)
    xx, yy = np.meshgrid(coordinates, coordinates)
    jones = np.broadcast_to(np.eye(2), (3, 3, 2, 2)).copy()
    mask = np.hypot(xx, yy) <= 1
    axis, artists = plot_polarization_ellipses_across_field(
        xx, yy, jones, [1, 0], mask
    )
    assert len(artists) == int(mask.sum())
    plt.close(axis.figure)
