# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

import numpy as np

from prtlab import example_cell_phone_lens_ex9_system, polarization_ray_trace


def test_cell_phone_lens_forward_trace_matches_matlab():
    output = polarization_ray_trace(
        example_cell_phone_lens_ex9_system(),
        [0, 0, 1],
        [0, 0.741676704231623, -0.1],
        [0, 1],
    )
    assert len(output.interactions) == 11
    assert len(output.final_ray_ids) == 1
    expected_positions = np.array(
        [
            [0, 0.741676704231623, 0.159396163650448],
            [0, 0.645230541985311, 0.877580313238837],
            [0, 0.494627159228190, 1.514145877780190],
            [0, 0.475744268443730, 2.099866961680170],
            [0, 0.462203382394511, 2.206523032508730],
            [0, 0.390695188763493, 2.902540866859470],
            [0, 0.370728806122263, 3.038099081352560],
            [0, 0.224884251196325, 3.859665156164050],
            [0, 0.101756289354252, 4.540000000000000],
            [0, 0.0662688237287398, 4.840000000000000],
            [0, -0.00249861353302133, 5.219969610400000],
        ]
    )
    positions = np.array(
        [interaction.position.real for interaction in output.interactions]
    )
    np.testing.assert_allclose(positions, expected_positions, atol=1e-9)

    final_ray = output.rays[output.final_ray_ids[0] - 1]
    np.testing.assert_allclose(
        final_ray.k.real,
        [0, -0.178088336099878, 0.984014504235165],
        atol=1e-9,
    )
    np.testing.assert_allclose(
        final_ray.flux, 0.605123991795796, atol=1e-10
    )
