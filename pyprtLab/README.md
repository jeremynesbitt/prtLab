# prtLab Python

prtLab is a polarization ray-tracing toolkit for optical systems containing
isotropic and uniaxial materials. This directory contains the installable
Python implementation. The MATLAB implementation in the repository root
remains the reference and the python is a port of that implementation that preserves numerical parity.

## Installation

From the `pyprtLab` directory, install the package into an activated Python
environment:

```bash
python -m pip install .
```

To work on the package itself, use an editable installation with the
development tools:

```bash
python -m pip install -e ".[dev]"
pytest
```

The package name and import name differ only in capitalization. Import it as:

```python
import prtlab
```

## Quick Start

This traces the quarter-wave plate used in Chipman Example 20.1 and coherently
recombines its ordinary and extraordinary branches:

```python
from prtlab import (
    coherent_final_jones,
    jones_retardance_waves,
    polarization_ray_trace,
    quarter_wave_plate_system,
)

system = quarter_wave_plate_system()
result = polarization_ray_trace(
    system,
    [0, 0, 1],
    [0, 0, 0],
    [0, 1],
    {"encodePropagationPhaseInP": True},
)
jones = coherent_final_jones(result)
print("Retardance [waves]:", jones_retardance_waves(jones))
```

## Examples and Notebooks

The examples are the best entry point for learning the current API. Python
scripts are stored in [`examples`](examples), and direct translations of the
MATLAB notebooks are stored in [`docs/source`](docs/source):

- [`QuarterWavePlateAnalysis.ipynb`](docs/source/QuarterWavePlateAnalysis.ipynb)
- [`CellPhoneLensSystemExample.ipynb`](docs/source/CellPhoneLensSystemExample.ipynb)
- [`CassegrainTelescopeExample.ipynb`](docs/source/CassegrainTelescopeExample.ipynb)
- [`FresnelRhombAnalysis.ipynb`](docs/source/FresnelRhombAnalysis.ipynb)

For interactive notebook use, register the project environment once:

```bash
python -m ipykernel install --user \
    --name prtlab-python --display-name "prtLab Python"
```

Then select **prtLab Python** from Jupyter's kernel menu.

## Solid Scene Example

```python
import matplotlib.pyplot as plt
import numpy as np

from prtlab import (
    create_fresnel_rhomb_solid,
    plot_prt_system_3d,
    polarization_ray_trace_solid,
)

rhomb = create_fresnel_rhomb_solid(
    1.495,
    length=14,
    height=8,
    width=8,
    shear=11.02085,
    wavelength=0.589,
)
result = polarization_ray_trace_solid(
    rhomb,
    [0, np.sin(np.deg2rad(0.1)), np.cos(np.deg2rad(0.1))],
    [0, 5, -2],
    [1, 1] / np.sqrt(2),
    {"maxInteractions": 8, "minAmplitude": 1e-4},
)

plot_prt_system_3d(
    rhomb,
    result,
    {"PostExtend": 4, "ShowLegend": True},
)
plt.show()
```

The 3D renderer colors ordinary, extraordinary, and reflected branches
separately. Red glyphs show the real time-domain locus of each complex
electric field.

## Validation

Run the Python tests and notebook checks with:

```bash
pytest
python tools/validate_notebooks.py
python tools/validate_notebooks.py --execute
python tools/smoke_test_wheel.py
```

See [`PORTING_STATUS.md`](PORTING_STATUS.md) for current MATLAB-to-Python
coverage.

## MATLAB Synchronization

Every MATLAB function in `prtLab/` is classified in
`tools/matlab_port_manifest.json` as parity-tested, internal, deferred, or
legacy. Validate the inventory and show deferred APIs with:

```bash
python tools/audit_api_parity.py --list-deferred
```

After changing MATLAB code, report files changed since the last reviewed
Python synchronization point with:

```bash
python tools/check_matlab_sync.py
```

After those changes have been reviewed, ported, and tested, record the current
commit with:

```bash
python tools/check_matlab_sync.py --mark-current
```

The marker is a review checkpoint, not an automatic claim of numerical parity.

## License

prtLab is distributed under the BSD 3-Clause License. See [`LICENSE`](LICENSE).
