# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Trace one meridional ray through Chipman's cell-phone lens example 3."""

from prtlab import example_cell_phone_lens_ex3_system, polarization_ray_trace


def run_example():
    system = example_cell_phone_lens_ex3_system()
    return polarization_ray_trace(
        system, [0, 0, 1], [0, 0.4, -0.1], [0, 1],
        {"minAmplitude": 0.2},
    )


if __name__ == "__main__":
    output = run_example()
    final = output.rays[output.final_ray_ids[0] - 1]
    print("Interactions:", len(output.interactions))
    print("Final position:", final.position.real)
    print("Final direction:", final.k.real)
    print("Final flux:", final.flux)
