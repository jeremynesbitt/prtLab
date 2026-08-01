# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Parity and smoke tests for the remaining catalog examples."""

import numpy as np
from parity import load_json_baseline

from prtlab import (
    create_right_triangular_prism_solid,
    example_cell_phone_lens_ex3_system,
    example_plane_mirror_system,
    plot_prt_system_3d,
    polarization_ray_trace,
)

BASELINE = load_json_baseline("remainingExampleBaselines.json")
TOLERANCE = BASELINE["absoluteTolerance"]


def unpack(value):
    return np.asarray(value["real"]) + 1j * np.asarray(value["imag"])


def test_cell_phone_example_3_matches_matlab():
    result = polarization_ray_trace(
        example_cell_phone_lens_ex3_system(),
        [0, 0, 1], [0, 0.4, -0.1], [0, 1],
        {"minAmplitude": 0.2},
    )
    reference = BASELINE["cellPhoneExample3"]
    assert len(result.final_ray_ids) == 1
    final = result.rays[result.final_ray_ids[0] - 1]
    dominant_positions = []
    for ray_id in final.metadata["history"]:
        interaction = next(
            (
                item
                for item in result.interactions
                if item.incident.id == ray_id
            ),
            None,
        )
        if interaction is not None:
            dominant_positions.append(interaction.position)
    np.testing.assert_allclose(
        dominant_positions,
        unpack(reference["dominantPositions"]),
        atol=TOLERANCE,
    )
    np.testing.assert_allclose(final.k, unpack(reference["finalK"]), atol=TOLERANCE)
    np.testing.assert_allclose(final.field_e, unpack(reference["finalFieldE"]), atol=TOLERANCE)
    np.testing.assert_allclose(final.P, unpack(reference["finalP"]), atol=TOLERANCE)
    np.testing.assert_allclose(final.flux, unpack(reference["finalFlux"]), atol=TOLERANCE)
    np.testing.assert_allclose(final.opl, unpack(reference["finalOPL"]), atol=TOLERANCE)


def test_plane_mirror_matches_matlab_complex_fresnel_solution():
    reference = BASELINE["planeMirror"]
    angle = np.deg2rad(reference["angleDeg"])
    result = polarization_ray_trace(
        example_plane_mirror_system(),
        [0, np.sin(angle), np.cos(angle)], [0, 0, 0],
        np.array([1, 1]) / np.sqrt(2),
    )
    interaction = result.interactions[0]
    reflected = next(
        child for child in interaction.children
        if child.branch_type == "reflected"
    )
    np.testing.assert_allclose(reflected.k, unpack(reference["reflectedK"]), atol=TOLERANCE)
    np.testing.assert_allclose(reflected.field_e, unpack(reference["reflectedFieldE"]), atol=TOLERANCE)
    np.testing.assert_allclose(reflected.field_h, unpack(reference["reflectedFieldH"]), atol=TOLERANCE)
    np.testing.assert_allclose(reflected.P, unpack(reference["reflectedP"]), atol=TOLERANCE)
    np.testing.assert_allclose(reflected.flux, unpack(reference["reflectedFlux"]), atol=TOLERANCE)
    np.testing.assert_allclose(
        interaction.coefficients["As"][2],
        unpack(reference["As"])[2],
        atol=TOLERANCE,
    )
    np.testing.assert_allclose(
        interaction.coefficients["Ap"][3],
        unpack(reference["Ap"])[3],
        atol=TOLERANCE,
    )


def test_right_triangular_crystal_prism_geometry_only_plot():
    prism = create_right_triangular_prism_solid(
        {"nO": 1.6557, "nE": 1.4852},
        leg_y=12, leg_z=12, width=8, optic_axis=[1, 0, 0],
        wavelength=0.633, wavelength_units="um",
        name="Calcite right triangular prism",
    )
    assert [face.name for face in prism.faces] == [
        "z leg", "hypotenuse", "y leg"
    ]
    np.testing.assert_allclose(
        prism.material.axis_data["opticAxis"], [1, 0, 0], atol=1e-14
    )
    axis, handles = plot_prt_system_3d(
        prism, None, {"ShowRays": False, "ShowPolarization": False}
    )
    assert handles["surfaces"]
    assert not handles["rays"]
    axis.figure.clf()
