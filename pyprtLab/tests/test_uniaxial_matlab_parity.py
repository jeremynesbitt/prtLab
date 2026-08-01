# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import csv
from pathlib import Path

import numpy as np
import pytest
from parity import load_json_baseline

from prtlab import add_surface, create_optical_system, polarization_ray_trace
from prtlab.linalg import calc_k_from_theta_x_theta_y

BASELINE = Path(__file__).parent / "baselines" / "uniaxialRayTraceBaseline.csv"
N_O = 1.666450305378000
N_E = 1.490046084791000
METADATA = load_json_baseline("uniaxialRayTraceBaseline.meta.json")


def baseline_rows():
    with BASELINE.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def interface_system(optic_axis):
    system = create_optical_system(0.633)
    add_surface(
        system, np.inf, 0, "plane", {}, "isotropic", {"n": 1.0}, {},
        {"type": "bare"},
    )
    add_surface(
        system, np.inf, 0, "plane", {}, "uniaxial",
        {"nO": N_O, "nE": N_E}, {"opticAxis": optic_axis},
        {"type": "bare"},
    )
    return system


@pytest.mark.parametrize("row", baseline_rows(), ids=lambda row: row["caseId"])
def test_isotropic_to_uniaxial_matches_matlab_baseline(row):
    axis = np.array(
        [float(row["axisX"]), float(row["axisY"]), float(row["axisZ"])]
    )
    axis /= np.linalg.norm(axis)
    k_in = calc_k_from_theta_x_theta_y(
        float(row["thetaX_deg"]), float(row["thetaY_deg"])
    )
    result = polarization_ray_trace(
        interface_system(axis), k_in, [0, 0, 0], [1, 0]
    )
    interaction = result.interactions[0]
    children = {child.mode: child for child in interaction.children}
    ordinary = children["ordinary"]
    extraordinary = children["extraordinary"]

    expected_k_o = values(row, "kOrd")
    expected_k_e = values(row, "kExt")
    expected_s_o = values(row, "SOrd")
    expected_s_e = values(row, "SExt")
    tolerance = METADATA["absoluteTolerance"]
    np.testing.assert_allclose(ordinary.k.real, expected_k_o, atol=tolerance)
    np.testing.assert_allclose(extraordinary.k.real, expected_k_e, atol=tolerance)
    np.testing.assert_allclose(ordinary.S.real, expected_s_o, atol=tolerance)
    np.testing.assert_allclose(extraordinary.S.real, expected_s_e, atol=tolerance)

    np.testing.assert_allclose(ordinary.metadata["n"].real, float(row["nOrd"]),
                               atol=1e-12)
    np.testing.assert_allclose(extraordinary.metadata["n"].real,
                               float(row["nExt"]), atol=1e-12)
    np.testing.assert_allclose(extraordinary.metadata["n_SE"].real,
                               float(row["nExt_SE"]), atol=1e-12)
    np.testing.assert_allclose(ordinary.amplitude, float(row["ampOrd"]),
                               atol=1e-12)
    np.testing.assert_allclose(extraordinary.amplitude, float(row["ampExt"]),
                               atol=1e-12)
    np.testing.assert_allclose(ordinary.flux, float(row["fluxOrd"]),
                               atol=1e-12)
    np.testing.assert_allclose(extraordinary.flux, float(row["fluxExt"]),
                               atol=1e-12)

    # Transmitted mode coefficients agree up to eigenvector sign.
    for name in ("a_s_to", "a_s_te", "a_p_to", "a_p_te"):
        np.testing.assert_allclose(
            abs(interaction.coefficients[name]), abs(float(row[name])),
            atol=1e-12,
        )

    # MATLAB obtains its reflected pair from an SVD, so its two basis vectors
    # can be rotated relative to deterministic physical s/p. Column norms are
    # invariant to that output-basis rotation.
    for incident_mode in ("s", "p"):
        actual = np.hypot(
            abs(interaction.coefficients[f"a_{incident_mode}_rs"]),
            abs(interaction.coefficients[f"a_{incident_mode}_rp"]),
        )
        expected = np.hypot(
            float(row[f"a_{incident_mode}_rs"]),
            float(row[f"a_{incident_mode}_rp"]),
        )
        np.testing.assert_allclose(actual, expected, atol=1e-12)
    np.testing.assert_allclose(
        interaction.diagnostics["boundaryResidual_s"], 0, atol=1e-12
    )
    np.testing.assert_allclose(
        interaction.diagnostics["boundaryResidual_p"], 0, atol=1e-12
    )


def values(row, prefix):
    return np.array([float(row[prefix + component]) for component in "XYZ"])
