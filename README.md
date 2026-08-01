# prtLab

prtLab is a MATLAB-primary toolkit for polarization ray tracing of optical
systems such as waveplates, crystal polarizers, Fresnel rhombs, and reflective
telescopes. It is not a commercial optical-design program, but it is usable as
a research and educational tool.

The project was inspired by *Polarized Light and Optical Systems* by Chipman,
Lam, and Young (CRC Press, ISBN 9781498700566). It implements the polarization
ray-tracing framework developed in that text, including three-dimensional
electric-field propagation and polarization ray-tracing matrices.

The best way to engage with prtLab is to start with an example and build from
it. The examples exercise the same public interfaces used by the tests and the
longer notebook-based documentation.

## Installation

### MATLAB

Clone or download the repository, then add the MATLAB source and examples to
your path:

```matlab
addpath('prtLab');
addpath('examples');
```

Most ray tracing uses base MATLAB. The achromatic-waveplate thickness
optimization example requires Optimization Toolbox because it uses
`lsqnonlin`.

### Python

The Python port is distributed from the `pyprtLab` directory. From the
repository root, install it into an activated Python environment with:

```bash
python -m pip install ./pyprtLab
```

The import name is lowercase:

```python
import prtlab
```

MATLAB remains the reference implementation while numerical parity of the
Python port is established.

## Start With Examples

| Topic | MATLAB example | Detailed documentation | Python notebook |
| --- | --- | --- | --- |
| Quarter-wave plate | [`testQuarterWavePlate.m`](examples/testQuarterWavePlate.m) | [`QuarterWavePlateAnalysis.md`](docs/QuarterWavePlateAnalysis.md) | [`QuarterWavePlateAnalysis.ipynb`](pyprtLab/docs/source/QuarterWavePlateAnalysis.ipynb) |
| Cell-phone lens | [`testCellPhoneSystem.m`](examples/testCellPhoneSystem.m) | [`CellPhoneLensSystemExample.md`](docs/CellPhoneLensSystemExample.md) | [`CellPhoneLensSystemExample.ipynb`](pyprtLab/docs/source/CellPhoneLensSystemExample.ipynb) |
| Cassegrain telescope | [`testCassegrainTelescope.m`](examples/testCassegrainTelescope.m) | [`CassegrainTelescopeExample.md`](docs/CassegrainTelescopeExample.md) | [`CassegrainTelescopeExample.ipynb`](pyprtLab/docs/source/CassegrainTelescopeExample.ipynb) |
| Fresnel rhomb | [`testFresnelRhomb.m`](examples/testFresnelRhomb.m) | [`FresnelRhombExample.md`](docs/FresnelRhombExample.md) | [`FresnelRhombAnalysis.ipynb`](pyprtLab/docs/source/FresnelRhombAnalysis.ipynb) |

Additional MATLAB examples cover achromatic and optically contacted
waveplates, Glan-Taylor polarizers, plane mirrors, and multi-solid scenes. The
Python scripts in [`pyprtLab/examples`](pyprtLab/examples) are direct ports of
the corresponding MATLAB examples where available.

## MATLAB Quick Start

This traces a quarter-wave plate as a bare bones test of the method. Propagation
phase is encoded in the accumulated polarization matrices so the two crystal
eigenmodes can be coherently recombined.

```matlab
T = exampleQuarterWavePlateSystem();
k = [0; 0; 1];
options = struct('encodePropagationPhaseInP', true);
rayOutput = polarizationRayTrace(T, k, [0; 0; 0], [0; 1], options);

Ptot = zeros(3, 3);
for rayId = rayOutput.finalRayIds
    Ptot = Ptot + rayOutput.rays(rayId).P;
end

Eout = Ptot * [0; 1; 0]
```

## Python Quick Start

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

## Current Scope

prtLab currently provides:

- isotropic and uniaxial materials;
- isotropic/isotropic, isotropic/uniaxial, uniaxial/isotropic, and
  uniaxial/uniaxial interfaces;
- planar, spherical, and odd-asphere sequential surfaces;
- extruded solids and multi-solid scenes;
- transmitted and reflected ray branches, including total internal reflection;
- optical-path-length and polarization-matrix accumulation;
- Jones analysis and two- and three-dimensional plotting helpers; and
- MATLAB/Python numerical-parity tests for core calculations and examples.

The built-in ray tracer is bare bones and is not a replacement
for a full optical-design ray tracer.  Coatings and biaxial materials are not yet
fully implemented.

## Project Structure

- `prtLab/`: MATLAB source
- `examples/`: MATLAB systems and worked examples
- `tests/`: MATLAB physics and regression tests
- `docs/`: generated Markdown and source MATLAB notebooks
- `pyprtLab/`: installable Python package, tests, examples, and notebooks

## Future Work

Planned areas include additional crystal-polarizer examples, biaxial-material
support, expanded coating models, and optional integration with a more complete
geometric ray tracer. Suggestions and reproducible test cases are welcome
through GitHub issues.

## License

prtLab is distributed under the BSD 3-Clause License. See [`LICENSE`](LICENSE).
