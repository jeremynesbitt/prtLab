# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import numpy as np

from prtlab import (
    achromatic_retardance_spectrum,
    mgf2_optical_constants,
    n_pk52a_optical_constants,
    quartz_optical_constants,
    scalar_crossed_axis_retardance,
)


def test_reference_dispersion_values_at_633_nm():
    np.testing.assert_allclose(
        quartz_optical_constants(0.633),
        [1.5425991960551284, 1.5516438611680567],
        atol=2e-15,
    )
    np.testing.assert_allclose(
        mgf2_optical_constants(0.633),
        [1.3769811265653429, 1.3887595196657687],
        atol=2e-15,
    )
    np.testing.assert_allclose(
        n_pk52a_optical_constants(0.633), 1.49570788594985,
        atol=2e-15,
    )


def test_vectorized_dispersion_preserves_shape():
    wavelengths = np.linspace(0.4, 0.7, 7)
    n_o, n_e = quartz_optical_constants(wavelengths)
    assert n_o.shape == wavelengths.shape
    assert n_e.shape == wavelengths.shape


def test_prt_and_scalar_achromatic_retardance_agree_at_normal_incidence():
    wavelengths = np.linspace(0.4, 0.7, 7)
    traced = achromatic_retardance_spectrum(
        wavelengths, 230.537086, 190.217463
    )
    scalar = scalar_crossed_axis_retardance(
        wavelengths, 230.537086, 190.217463
    )
    np.testing.assert_allclose(traced, scalar, atol=2e-12)
