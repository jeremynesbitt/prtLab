# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Trace symmetric marginal rays through the Cassegrain telescope."""

from prtlab import (
    example_cassegrain_telescope,
    plot_prt_system_3d,
    polarization_ray_trace,
)


def run_example():
    system = example_cassegrain_telescope()
    upper = polarization_ray_trace(system, [0, 0, 1], [0, 4000, -6000], [1, 0])
    lower = polarization_ray_trace(system, [0, 0, 1], [0, -4000, -6000], [1, 0])
    return system, [upper, lower]


def plot_example():
    system, outputs = run_example()
    return plot_prt_system_3d(
        system, outputs,
        {"ShowPolarization": False, "EqualAxes": False, "PostExtend": 500},
    )


if __name__ == "__main__":
    _, outputs = run_example()
    for result in outputs:
        ray = result.rays[result.final_ray_ids[0] - 1]
        print("Image position:", ray.position.real, "OPL:", ray.opl)
