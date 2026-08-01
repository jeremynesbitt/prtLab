# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Polarization ray tracing for simple optical systems."""

from ._version import __version__
from .analysis import (
    coherent_final_jones,
    coherent_final_p,
    jones_major_axis_ellipticity,
    jones_retardance,
    jones_retardance_waves,
    retardance_by_qp_method,
    select_dominant_final_ray,
)
from .catalog import (
    example_cassegrain_telescope,
    example_cell_phone_lens_ex3_system,
    example_cell_phone_lens_ex9_system,
    example_plane_mirror_system,
)
from .coordinates import (
    calc_o,
    dipole_basis_vectors,
    double_pole_basis_vectors,
    sp_basis,
    transform_p_to_jones,
)
from .dispersion import (
    mgf2_optical_constants,
    n_pk52a_optical_constants,
    quartz_optical_constants,
)
from .linalg import (
    calc_k_from_theta_x_theta_y,
    get_theta_x_theta_y_from_k,
    make_k,
    normalize,
    prt_norm,
    snell_vector,
)
from .models import (
    Interaction,
    Medium,
    OpticalSystem,
    RayBranch,
    Scene,
    Solid,
    SolidFace,
    Surface,
    TraceOptions,
    TraceResult,
)
from .plotting import (
    PlotOptions,
    plot_basis_vectors_on_sphere,
    plot_jones_pupil,
    plot_jones_vector,
    plot_polarization_ellipses_across_field,
    plot_polarization_glyph_3d,
    plot_prt_lens_cross_section,
    plot_prt_ray_trace,
    plot_prt_system_3d,
    polarization_ellipse_points_2d,
    prepare_jones_pupil,
)
from .scene_tracing import polarization_ray_trace_scene, polarization_ray_trace_solid
from .scenes import add_solid, create_glan_taylor_polarizer_scene, create_scene
from .solids import (
    create_extruded_polygon_solid,
    create_fresnel_rhomb_solid,
    create_right_triangular_prism_solid,
    isotropic_medium,
    material_medium,
    transform_solid,
)
from .systems import (
    add_surface,
    create_optical_system,
    update_clear_apertures_from_ray_trace,
)
from .tracing import polarization_ray_trace
from .waveplates import (
    achromatic_retardance_spectrum,
    achromatic_wave_plate_system,
    optimize_achromatic_wave_plate,
    quarter_wave_plate_system,
    scalar_crossed_axis_retardance,
)

__all__ = [
    "Interaction",
    "Medium",
    "OpticalSystem",
    "PlotOptions",
    "RayBranch",
    "Scene",
    "Solid",
    "SolidFace",
    "Surface",
    "TraceOptions",
    "TraceResult",
    "__version__",
    "achromatic_retardance_spectrum",
    "achromatic_wave_plate_system",
    "add_solid",
    "add_surface",
    "calc_k_from_theta_x_theta_y",
    "calc_o",
    "coherent_final_jones",
    "coherent_final_p",
    "create_extruded_polygon_solid",
    "create_fresnel_rhomb_solid",
    "create_glan_taylor_polarizer_scene",
    "create_optical_system",
    "create_right_triangular_prism_solid",
    "create_scene",
    "dipole_basis_vectors",
    "double_pole_basis_vectors",
    "example_cassegrain_telescope",
    "example_cell_phone_lens_ex3_system",
    "example_cell_phone_lens_ex9_system",
    "example_plane_mirror_system",
    "get_theta_x_theta_y_from_k",
    "isotropic_medium",
    "jones_major_axis_ellipticity",
    "jones_retardance",
    "jones_retardance_waves",
    "make_k",
    "material_medium",
    "mgf2_optical_constants",
    "n_pk52a_optical_constants",
    "normalize",
    "optimize_achromatic_wave_plate",
    "plot_basis_vectors_on_sphere",
    "plot_jones_pupil",
    "plot_jones_vector",
    "plot_polarization_ellipses_across_field",
    "plot_polarization_glyph_3d",
    "plot_prt_lens_cross_section",
    "plot_prt_ray_trace",
    "plot_prt_system_3d",
    "polarization_ellipse_points_2d",
    "polarization_ray_trace",
    "polarization_ray_trace_scene",
    "polarization_ray_trace_solid",
    "prepare_jones_pupil",
    "prt_norm",
    "quarter_wave_plate_system",
    "quartz_optical_constants",
    "retardance_by_qp_method",
    "scalar_crossed_axis_retardance",
    "select_dominant_final_ray",
    "snell_vector",
    "sp_basis",
    "transform_p_to_jones",
    "transform_solid",
    "update_clear_apertures_from_ray_trace",
]
