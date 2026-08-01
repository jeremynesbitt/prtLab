#!/usr/bin/env python3
# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Build, install, and exercise a wheel outside the source tree."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(*arguments: str, cwd: Path = ROOT, env=None) -> None:
    subprocess.run(arguments, cwd=cwd, env=env, check=True)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="prtlab-wheel-") as temporary:
        temporary_path = Path(temporary)
        wheel_directory = temporary_path / "wheel"
        install_directory = temporary_path / "installed"
        run(
            sys.executable, "-m", "pip", "wheel", ".",
            "--no-deps", "--no-build-isolation", "--wheel-dir",
            str(wheel_directory),
        )
        wheels = list(wheel_directory.glob("*.whl"))
        if len(wheels) != 1:
            raise RuntimeError(f"expected one wheel, found {wheels}")
        run(
            sys.executable, "-m", "pip", "install", str(wheels[0]),
            "--no-deps", "--target", str(install_directory),
        )
        environment = os.environ.copy()
        environment["PYTHONPATH"] = str(install_directory)
        smoke_code = """
import prtlab
from prtlab import polarization_ray_trace, quarter_wave_plate_system

assert prtlab.__version__ == "0.1.0.dev0"
result = polarization_ray_trace(
    quarter_wave_plate_system(), [0, 0, 1], [0, 0, 0], [0, 1]
)
assert len(result.final_ray_ids) == 2
print(f"Installed prtLab {prtlab.__version__}; final rays: {result.final_ray_ids}")
"""
        run(
            sys.executable, "-c", smoke_code,
            cwd=temporary_path, env=environment,
        )


if __name__ == "__main__":
    main()
