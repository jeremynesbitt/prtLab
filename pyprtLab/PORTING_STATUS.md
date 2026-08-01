# Porting Status

The MATLAB implementation is authoritative until this file marks a subsystem
as parity-tested.

`tools/matlab_port_manifest.json` is the machine-checked function inventory.
All currently classified public MATLAB helpers have Python equivalents.

| MATLAB subsystem | Python subsystem | Status |
| --- | --- | --- |
| `createOpticalSystem.m`, `addSurface.m` | `prtlab.models` | Implemented |
| Vector and coordinate helpers | `prtlab.linalg`, `prtlab.coordinates` | P-to-Jones coordinate systems MATLAB parity-tested |
| Jones ellipse analysis and pupil plotting | `prtlab.analysis`, `prtlab.plotting` | Implemented; ellipse analysis MATLAB parity-tested |
| `traceIsotropicToIsotropic.m` | `prtlab.interfaces.isotropic` | Implemented |
| `traceIsotropicToUniaxial.m` | `prtlab.interfaces.uniaxial` | MATLAB parity-tested |
| `polarizationRayTrace.m` | `prtlab.tracing` | Planar sequential slice |
| `traceUniaxialToIsotropic.m` | `prtlab.interfaces.uniaxial` | MATLAB parity-tested |
| `traceUniaxialToUniaxial.m` | `prtlab.interfaces.uniaxial` | MATLAB parity-tested |
| Plane and odd-asphere intersection | `prtlab.geometry` | Implemented |
| Extruded polygon solids and rigid transforms | `prtlab.solids` | Implemented |
| `polarizationRayTraceScene.m` | `prtlab.scene_tracing` | MATLAB parity-tested for rhomb and Glan-Taylor scenes |
| Fresnel rhomb and two-prism examples | `examples/fresnel_rhomb.py`, `examples/two_prism_scene.py` | Fresnel rhomb MATLAB parity-tested |
| Glan-Taylor polarizer scene | `examples/glan_taylor_polarizer.py` | MATLAB parity-tested |
| `plotPrtSystem3D.m` | `prtlab.plotting.plot_prt_system_3d` | Initial solid/scene port |
| Quarter-wave plate example | `examples/quarter_wave_plate.py` | MATLAB parity-tested |
| Contacted waveplate example | `examples/contacted_wave_plate.py` | MATLAB parity-tested |
| Quartz, MgF2, and N-PK52A dispersion | `prtlab.dispersion` | Implemented |
| Achromatic waveplate optimization | `prtlab.waveplates` | Implemented |
| Cell-phone lens example 9 | `examples/cell_phone_lens.py` | MATLAB parity-tested |
| Cell-phone lens example 3 | `examples/cell_phone_lens_ex3.py` | Dominant-path MATLAB parity-tested |
| Plane mirror example | `examples/plane_mirror.py` | MATLAB parity-tested |
| Right triangular crystal prism | `examples/right_triangular_crystal_prism.py` | Implemented and geometry-tested |
| Quarter-waveplate Python notebook | `docs/source/QuarterWavePlateAnalysis.ipynb` | Implemented |
| Achromatic-waveplate Python notebook | `docs/source/AchromaticWavePlateAnalysis.ipynb` | Implemented |
| Cell-phone lens Python notebook | `docs/source/CellPhoneLensSystemExample.ipynb` | Implemented |
| Reflected sequential routing | `prtlab.tracing` | MATLAB Cassegrain parity-tested |
| Cassegrain Python notebook | `docs/source/CassegrainTelescopeExample.ipynb` | Implemented |
| Fresnel-rhomb Python notebook | `docs/source/FresnelRhombAnalysis.ipynb` | Implemented and MATLAB parity-tested |
| Glan-Taylor Python notebook | `docs/source/GlanTaylorPolarizerAnalysis.ipynb` | Implemented and MATLAB parity-tested |

Numerical parity fixtures cover isotropic-to-uniaxial ray grids and detailed
normal/oblique uniaxial-to-isotropic and uniaxial-to-uniaxial interactions.
Scene fixtures cover repeated-TIR Fresnel-rhomb paths and every surviving
Glan-Taylor branch, including accumulated fields, P/Q matrices, flux, and OPL.

MATLAB baseline provenance is recorded beside each fixture. Run
`tools/check_matlab_sync.py` to detect MATLAB changes after the reviewed sync
commit. The GitHub Actions Python matrix validates lint, tests, notebooks, and
an installed wheel on supported Python versions.
