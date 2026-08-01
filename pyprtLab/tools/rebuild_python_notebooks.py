#!/usr/bin/env python3
# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Rebuild Python notebooks as faithful translations of MATLAB notebooks."""

from __future__ import annotations

import copy
import json
import textwrap
from pathlib import Path

PYTHON_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PYTHON_ROOT.parent
MATLAB_NOTEBOOKS = REPOSITORY_ROOT / "docs" / "source"
PYTHON_NOTEBOOKS = PYTHON_ROOT / "docs" / "source"
KERNELSPEC = {
    "display_name": "prtLab Python",
    "language": "python",
    "name": "prtlab-python",
}


def source(code: str) -> list[str]:
    text = textwrap.dedent(code).strip("\n") + "\n"
    return text.splitlines(keepends=True)


TRANSLATIONS = {
    "QuarterWavePlateAnalysis.ipynb": (
        "QuarterWavePlateAnalysis.ipynb",
        [
            """
            import copy

            import matplotlib.pyplot as plt
            import numpy as np

            from prtlab import (
                add_surface,
                calc_k_from_theta_x_theta_y,
                create_optical_system,
                plot_polarization_ellipses_across_field,
                plot_prt_lens_cross_section,
                plot_prt_ray_trace,
                polarization_ray_trace,
                transform_p_to_jones,
                update_clear_apertures_from_ray_trace,
            )
            """,
            """
            wavelength = 0.633  # um
            n_o = 1.656
            n_e = 1.485
            t_qwp = 1.25 * wavelength / (n_o - n_e)
            optic_axis = np.array([-np.sqrt(2) / 2, -np.sqrt(2) / 2, 0])

            T = create_optical_system(wavelength)
            T.wavelength_units = "um"

            air_index = {"n": 1.0}
            calcite_index = {"nO": n_o, "nE": n_e}
            calcite_axis = {"opticAxis": optic_axis}
            bare_coating = {"type": "bare"}

            # Object/input space. The zero-thickness row makes the incident
            # medium explicit without adding propagation.
            add_surface(
                T, np.inf, 0, "plane", {}, "isotropic",
                air_index, {}, bare_coating,
            )

            # Plane entrance surface followed by propagation through the
            # calcite plate.
            add_surface(
                T, np.inf, t_qwp, "plane", {}, "uniaxial",
                calcite_index, calcite_axis, bare_coating,
            )

            # Plane exit surface into air. Thickness is zero because this row
            # describes the output medium after the waveplate.
            add_surface(
                T, np.inf, 0, "plane", {}, "isotropic",
                air_index, {}, bare_coating,
            )
            """,
            """
            T
            """,
            """
            k_in = np.array([0, 0, 1])  # normal incidence
            pos_in = np.array([0, 0, 0])  # on axis
            E_in = np.array([0, 1])  # can be 2D or 3D vector
            trace_options = {"encodePropagationPhaseInP": True}
            ray_output = polarization_ray_trace(
                T, k_in, pos_in, E_in, trace_options
            )
            """,
            """
            ray_output
            """,
            """
            first, second = ray_output.final_ray_ids[:2]
            opl = (
                ray_output.rays[first - 1].opl
                - ray_output.rays[second - 1].opl
            )
            opl_in_waves = opl / wavelength
            opl_in_waves
            """,
            """
            angles = np.linspace(0, 5, 21)
            opl_vs_aoi = np.zeros_like(angles)

            for index, angle in enumerate(angles):
                theta = np.deg2rad(angle)
                k_in = np.array([0, np.sin(theta), np.cos(theta)])
                ray_output = polarization_ray_trace(
                    T, k_in, pos_in, E_in, trace_options
                )
                first, second = ray_output.final_ray_ids[:2]
                opl = (
                    ray_output.rays[first - 1].opl
                    - ray_output.rays[second - 1].opl
                )
                opl_vs_aoi[index] = opl / wavelength

            plt.figure()
            _ = plt.plot(angles, opl_vs_aoi - 1)
            plt.xlabel("Angle of Incidence [deg]")
            _ = plt.ylabel("Retardance [waves]")
            """,
            """
            Tm = copy.deepcopy(T)
            # Integer multiple guarantees quarter wave at 0 AOI.
            Tm.surfaces[1].thickness = 17 * T.surfaces[1].thickness
            opl_multiwave_vs_aoi = np.zeros_like(angles)

            for index, angle in enumerate(angles):
                theta = np.deg2rad(angle)
                k_in = np.array([0, np.sin(theta), np.cos(theta)])
                ray_output = polarization_ray_trace(
                    Tm, k_in, pos_in, E_in, trace_options
                )
                first, second = ray_output.final_ray_ids[:2]
                opl = (
                    ray_output.rays[first - 1].opl
                    - ray_output.rays[second - 1].opl
                )
                opl_multiwave_vs_aoi[index] = opl / wavelength

            plt.figure()
            _ = plt.plot(angles, opl_multiwave_vs_aoi - 21)
            plt.xlabel("Angle of Incidence [deg]")
            _ = plt.ylabel("Retardance [waves]")
            """,
            """
            X, Y = np.meshgrid(np.linspace(-5, 5, 17), np.linspace(-5, 5, 17))
            tX, tY = np.meshgrid(np.linspace(-5, 5, 17), np.linspace(-5, 5, 17))
            coordinate = {
                "type": "doublePole",
                "a_loc": np.array([0, 0, 1]),
                "x_o": np.array([1, 0, 0]),
            }
            Jall = np.zeros((*X.shape, 2, 2), dtype=complex)
            Jall_m = np.zeros_like(Jall)

            for index in np.ndindex(X.shape):
                k_in = calc_k_from_theta_x_theta_y(tX[index], tY[index])
                pos_in = np.array([X[index], Y[index], 0])
                ray_output = polarization_ray_trace(
                    T, k_in, pos_in, E_in, trace_options
                )
                P_tot = sum(
                    (ray_output.rays[ray_id - 1].P
                     for ray_id in ray_output.final_ray_ids),
                    np.zeros((3, 3), dtype=complex),
                )
                k_out = ray_output.rays[ray_output.final_ray_ids[0] - 1].k
                Jall[index] = transform_p_to_jones(
                    P_tot, k_in, k_out, coordinate
                )

                ray_output_m = polarization_ray_trace(
                    Tm, k_in, pos_in, E_in, trace_options
                )
                P_tot_m = sum(
                    (ray_output_m.rays[ray_id - 1].P
                     for ray_id in ray_output_m.final_ray_ids),
                    np.zeros((3, 3), dtype=complex),
                )
                k_out_m = ray_output_m.rays[
                    ray_output_m.final_ray_ids[0] - 1
                ].k
                Jall_m[index] = transform_p_to_jones(
                    P_tot_m, k_in, k_out_m, coordinate
                )

            axis, _ = plot_polarization_ellipses_across_field(
                X, Y, Jall, E_in
            )
            _ = axis.set_title("True 0 Order Plate")
            axis, _ = plot_polarization_ellipses_across_field(
                X, Y, Jall_m, E_in
            )
            _ = axis.set_title("Multi Order plate")
            """,
            """
            tta = 70
            theta = np.deg2rad(tta)
            k_in = np.array([0, np.sin(theta), np.cos(theta)])
            pos_in = np.array([0, 0, 0])
            E_in = np.array([0, 1])
            ray_output = polarization_ray_trace(
                T, k_in, pos_in, E_in, trace_options
            )

            Tplot = copy.deepcopy(T)
            update_clear_apertures_from_ray_trace(
                Tplot, ray_output, margin=1.15, minimum=0.25
            )
            axis, _ = plot_prt_lens_cross_section(Tplot)
            _ = plot_prt_ray_trace(ray_output, ax=axis)
            """,
        ],
    ),
    "CellPhoneLensSystemExample.ipynb": (
        "CellPhoneLensSystemExample.ipynb",
        [
            """
            import copy

            import matplotlib.pyplot as plt
            import numpy as np

            from prtlab import (
                double_pole_basis_vectors,
                example_cell_phone_lens_ex3_system,
                plot_basis_vectors_on_sphere,
                plot_prt_lens_cross_section,
                plot_prt_ray_trace,
                polarization_ray_trace,
                sp_basis,
                transform_p_to_jones,
                update_clear_apertures_from_ray_trace,
            )
            """,
            """
            T = example_cell_phone_lens_ex3_system()
            T
            """,
            """
            Ein = np.array([1, 0])
            k = np.array([0, np.sin(np.deg2rad(30)), np.cos(np.deg2rad(30))])
            x1 = np.array([0, 0.775, 0])
            x2 = np.array([0, -0.775, 0])
            options = {"minRelativeFlux": 0.25}

            d1 = polarization_ray_trace(T, k, x1, Ein, options)
            d2 = polarization_ray_trace(T, k, x2, Ein, options)

            Tplot = copy.deepcopy(T)
            update_clear_apertures_from_ray_trace(
                Tplot, [d1, d2], margin=1.30, minimum=0.25
            )

            k2 = np.array([0, np.sin(np.deg2rad(20)), np.cos(np.deg2rad(20))])
            k3 = np.array([0, 0, 1])
            d3 = polarization_ray_trace(Tplot, k2, x1, Ein, options)
            d4 = polarization_ray_trace(Tplot, k2, x2, Ein, options)
            d5 = polarization_ray_trace(Tplot, k3, x1, Ein, options)
            d6 = polarization_ray_trace(Tplot, k3, x2, Ein, options)
            axis, _ = plot_prt_lens_cross_section(Tplot)
            _ = plot_prt_ray_trace([d1, d2, d3, d4, d5, d6], ax=axis)
            T = Tplot
            """,
            """
            pupil_axis = np.linspace(-0.775, 0.775, 51)
            X, Y = np.meshgrid(pupil_axis, pupil_axis)
            pupil = np.hypot(X, Y) <= 0.775
            E_in = np.array([1, 0])
            k_in = np.array([0, 0, 1])
            coordinate = {
                "type": "doublePole", "a_loc": k_in,
                "x_o": np.array([1, 0, 0]),
            }
            Jall = np.zeros((*X.shape, 2, 2), dtype=complex)
            I_out = np.zeros_like(X)

            for index in np.ndindex(X.shape):
                pos_in = np.array([X[index], Y[index], 0])
                ray_output = polarization_ray_trace(
                    T, k_in, pos_in, E_in, options
                )
                if ray_output.final_ray_ids:
                    final_ray = ray_output.rays[
                        ray_output.final_ray_ids[0] - 1
                    ]
                    J = transform_p_to_jones(
                        final_ray.P, k_in, final_ray.k, coordinate
                    )
                    Jall[index] = J
                    E_out = np.array([0, 1]) @ J @ E_in
                    I_out[index] = np.abs(E_out) ** 2

            plt.figure()
            plt.imshow(I_out * pupil)
            _ = plt.colorbar()
            """,
            """
            tta = 7
            k_in = np.array([0, np.sin(np.deg2rad(tta)), np.cos(np.deg2rad(tta))])
            pupil_axis = np.linspace(
                -0.741676704231623, 0.741676704231623, 17
            )
            X, Y = np.meshgrid(pupil_axis, pupil_axis)
            Ein = np.array([1, 0])
            pupil = np.hypot(X, Y) <= 0.741676704231623
            x_loc, y_loc, k_loc = [], [], []
            coordinate = {
                "type": "doublePole", "a_loc": k_in,
                "x_o": np.array([1, 0, 0]),
            }

            for index in np.ndindex(X.shape):
                ray_output = polarization_ray_trace(
                    T, k_in, [X[index], Y[index], 0], Ein, options
                )
                if not ray_output.final_ray_ids:
                    continue
                final_ray = ray_output.rays[ray_output.final_ray_ids[0] - 1]
                if (
                    np.max(np.abs(np.imag(final_ray.k))) < 1e-12
                    and pupil[index]
                ):
                    x_out, y_out = double_pole_basis_vectors(
                        final_ray.k, coordinate["a_loc"], coordinate["x_o"]
                    )
                    x_loc.append(x_out)
                    y_loc.append(y_out)
                    k_loc.append(final_ray.k)

            axis, _ = plot_basis_vectors_on_sphere(
                np.real(k_loc), np.real(x_loc), None,
                show_sphere=False, arrow_scale=0.02,
                line_width=1, view=(0, tta, -90),
            )
            axis.set_axis_off()
            """,
            """
            # This is the only difference from the last snippet.
            x_loc, y_loc, k_loc = [], [], []
            normal = np.array([0, 0, 1])

            for index in np.ndindex(X.shape):
                ray_output = polarization_ray_trace(
                    T, k_in, [X[index], Y[index], 0], Ein, options
                )
                if not ray_output.final_ray_ids:
                    continue
                final_ray = ray_output.rays[ray_output.final_ray_ids[0] - 1]
                if (
                    np.max(np.abs(np.imag(final_ray.k))) < 1e-12
                    and pupil[index]
                ):
                    p_out, s_out = sp_basis(final_ray.k, normal)
                    x_loc.append(s_out)
                    y_loc.append(p_out)
                    k_loc.append(final_ray.k)

            axis, _ = plot_basis_vectors_on_sphere(
                np.real(k_loc), np.real(x_loc), None,
                show_sphere=False, arrow_scale=0.02,
                line_width=1, view=(0, tta, -90),
            )
            axis.set_axis_off()
            """,
        ],
    ),
    "CassegrainTelescopeExample.ipynb": (
        "CassegrainTelescopeExample.ipynb",
        [
            """
            import copy

            import numpy as np

            from prtlab import (
                example_cassegrain_telescope,
                plot_jones_pupil,
                plot_prt_lens_cross_section,
                plot_prt_ray_trace,
                polarization_ray_trace,
                select_dominant_final_ray,
                transform_p_to_jones,
                update_clear_apertures_from_ray_trace,
            )
            """,
            """
            T = example_cassegrain_telescope()
            T
            """,
            """
            k0 = np.array([0, 0, 1])
            x0 = np.array([0, 4000, -6000])
            x1 = np.array([0, -4000, -6000])
            options = {"minAmplitude": 0.05}
            d1 = polarization_ray_trace(T, k0, x0, [1, 0], options)
            d2 = polarization_ray_trace(T, k0, x1, [1, 0], options)
            Tplot = copy.deepcopy(T)
            update_clear_apertures_from_ray_trace(Tplot, [d1, d2])
            axis, _ = plot_prt_lens_cross_section(Tplot)
            _ = plot_prt_ray_trace([d1, d2], ax=axis)
            """,
            """
            d = 4161.65
            pupil_axis = np.linspace(-d, d, 50)
            X, Y = np.meshgrid(pupil_axis, pupil_axis)
            rho = np.hypot(X, Y)
            ftr = (rho >= 1000) & (rho <= d)
            E_in = np.array([1, 0])
            k_in = np.array([0, 0, 1])
            Jall = np.full((*X.shape, 2, 2), np.nan + 0j)
            OPL = np.full(X.shape, np.nan)
            coordinate = {
                "type": "doublePole", "a_loc": k_in,
                "x_o": np.array([1, 0, 0]),
            }

            def trace_pupil(system):
                jones = np.full((*X.shape, 2, 2), np.nan + 0j)
                opl = np.full(X.shape, np.nan)
                for index in np.ndindex(X.shape):
                    ray_output = polarization_ray_trace(
                        system, k_in, [X[index], Y[index], -6000], E_in
                    )
                    if ray_output.final_ray_ids:
                        final_ray, _ = select_dominant_final_ray(ray_output)
                        jones[index] = transform_p_to_jones(
                            final_ray.P, k_in, final_ray.k, coordinate
                        )
                        opl[index] = final_ray.opl
                return jones, opl

            Jall, OPL = trace_pupil(T)
            _ = plot_jones_pupil(Jall, ftr)
            """,
            """
            T.surfaces[1].index_data["n"] = 1e15
            T.surfaces[2].index_data["n"] = 1e15
            Jall, _ = trace_pupil(T)
            _ = plot_jones_pupil(Jall, ftr)
            """,
            """
            n_gold = 0.8 + 1.8j
            T.surfaces[1].index_data["n"] = n_gold
            T.surfaces[2].index_data["n"] = n_gold
            Jall, _ = trace_pupil(T)
            _ = plot_jones_pupil(Jall, ftr)
            """,
        ],
    ),
    "FresnelRhombAnalysis.ipynb": (
        "FresnelRhombExample.ipynb",
        [
            """
            import matplotlib.pyplot as plt
            import numpy as np

            from prtlab import (
                create_fresnel_rhomb_solid,
                jones_retardance,
                n_pk52a_optical_constants,
                plot_prt_system_3d,
                polarization_ray_trace_solid,
                select_dominant_final_ray,
                transform_p_to_jones,
            )
            """,
            """
            def fresnel_tir_coefficients(theta_deg, n):
                theta = np.deg2rad(theta_deg)
                sqrt_term = np.sqrt(np.sin(theta) ** 2 - n**2)
                cosine = np.cos(theta)
                r_s = (cosine - 1j * sqrt_term) / (cosine + 1j * sqrt_term)
                r_p = (
                    (n**2 * cosine - 1j * sqrt_term)
                    / (n**2 * cosine + 1j * sqrt_term)
                )
                return r_s, r_p

            angles = np.linspace(43, 60, 21)
            phase_difference = np.zeros_like(angles)
            for index, angle in enumerate(angles):
                r_s, r_p = fresnel_tir_coefficients(angle, 1 / 1.5)
                phase_difference[index] = np.angle(r_s) - np.angle(r_p)

            plt.figure()
            _ = plt.plot(angles, phase_difference)
            """,
            """
            wavelength = 0.589  # um
            n_glass = 1.5
            rhomb = create_fresnel_rhomb_solid(
                n_glass, length=14, height=8, width=8, shear=11.33,
                wavelength=wavelength / 1000, wavelength_units="mm",
            )
            options = {"maxInteractions": 8, "minAmplitude": 1e-4}
            k_in = np.array([0, 0, 1])
            x_in = np.array([0, 5, -2])
            E_in = np.array([1, 1]) / np.sqrt(2)
            ray_output = polarization_ray_trace_solid(
                rhomb, k_in, x_in, E_in, options
            )
            axis, _ = plot_prt_system_3d(
                rhomb, ray_output,
                {"RaySelection": "dominant", "PostExtend": 4},
            )
            _ = axis.set_title("Fresnel rhomb polarization evolution")
            """,
            """
            shears = np.linspace(10, 12, 101)
            coordinate = {
                "type": "doublePole", "a_loc": np.array([0, 0, 1]),
                "x_o": np.array([1, 0, 0]),
            }
            wavelength = 0.589
            n_glass = n_pk52a_optical_constants(wavelength)
            retardance_vs_shear = np.zeros_like(shears)
            aoi = np.zeros_like(shears)

            for index, shear in enumerate(shears):
                rhomb = create_fresnel_rhomb_solid(
                    n_glass, length=14, height=8, width=8, shear=shear,
                    wavelength=wavelength, wavelength_units="um",
                )
                ray_output = polarization_ray_trace_solid(
                    rhomb, k_in, x_in, E_in, options
                )
                final_ray, _ = select_dominant_final_ray(ray_output)
                jones = transform_p_to_jones(
                    final_ray.P, k_in, final_ray.k, coordinate
                )
                retardance_vs_shear[index] = jones_retardance(jones)
                normal = ray_output.interactions[1].normal
                aoi[index] = np.rad2deg(
                    np.arccos(np.dot([0, 0, 1], normal).real)
                )

            plt.figure()
            _ = plt.plot(aoi, retardance_vs_shear)
            """,
            """
            shear = 10.84  # eyeball fit
            wavelengths = np.linspace(0.4, 0.7, 21)
            retardance_vs_wavelength = np.zeros_like(wavelengths)

            for index, wavelength in enumerate(wavelengths):
                n_glass = n_pk52a_optical_constants(wavelength)
                rhomb = create_fresnel_rhomb_solid(
                    n_glass, length=14, height=8, width=8, shear=shear,
                    wavelength=wavelength, wavelength_units="um",
                )
                ray_output = polarization_ray_trace_solid(
                    rhomb, k_in, x_in, E_in, options
                )
                final_ray, _ = select_dominant_final_ray(ray_output)
                jones = transform_p_to_jones(
                    final_ray.P, k_in, final_ray.k, coordinate
                )
                retardance_vs_wavelength[index] = jones_retardance(jones)

            plt.figure()
            _ = plt.plot(wavelengths, retardance_vs_wavelength)
            """,
        ],
    ),
}

SCRIPT_NOTEBOOKS = {}


MARKDOWN_REPLACEMENTS = {
    "This is implemented using the table type in matlab.": (
        "This is implemented using the OpticalSystem type in Python."
    ),
    "the trace*.m files": "the trace interface files",
}


def translated_markdown(text: str) -> str:
    for old, new in MARKDOWN_REPLACEMENTS.items():
        text = text.replace(old, new)
    return text


def rebuild(target_name: str, matlab_name: str, translated_code: list[str]) -> None:
    matlab_path = MATLAB_NOTEBOOKS / matlab_name
    target_path = PYTHON_NOTEBOOKS / target_name
    matlab = json.loads(matlab_path.read_text(encoding="utf-8"))
    previous = json.loads(target_path.read_text(encoding="utf-8"))
    code_index = 0
    cells = []
    for original in matlab["cells"]:
        cell = copy.deepcopy(original)
        if cell["cell_type"] == "markdown":
            text = translated_markdown("".join(cell["source"]))
            cell["source"] = text.splitlines(keepends=True)
        elif cell["cell_type"] == "code":
            cell["source"] = source(translated_code[code_index])
            cell["execution_count"] = None
            cell["outputs"] = []
            code_index += 1
        cells.append(cell)
    if code_index != len(translated_code):
        raise RuntimeError(
            f"{target_name}: used {code_index} of {len(translated_code)} code cells"
        )
    previous["cells"] = cells
    previous.setdefault("metadata", {})["kernelspec"] = KERNELSPEC
    previous["metadata"]["language_info"] = {"name": "python", "version": "3"}
    target_path.write_text(
        json.dumps(previous, indent=1, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


def rebuild_script_notebook(target_name: str, translated_cells) -> None:
    target_path = PYTHON_NOTEBOOKS / target_name
    previous = json.loads(target_path.read_text(encoding="utf-8"))
    cells = []
    for cell_type, content in translated_cells:
        if cell_type == "markdown":
            cells.append(
                {
                    "cell_type": "markdown",
                    "metadata": {},
                    "source": source(content),
                }
            )
        else:
            cells.append(
                {
                    "cell_type": "code",
                    "execution_count": None,
                    "metadata": {},
                    "outputs": [],
                    "source": source(content),
                }
            )
    previous["cells"] = cells
    previous.setdefault("metadata", {})["kernelspec"] = KERNELSPEC
    previous["metadata"]["language_info"] = {"name": "python", "version": "3"}
    target_path.write_text(
        json.dumps(previous, indent=1, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    for target_name, (matlab_name, translated_code) in TRANSLATIONS.items():
        rebuild(target_name, matlab_name, translated_code)
        print(f"Rebuilt {target_name} from {matlab_name}")
    for target_name, translated_cells in SCRIPT_NOTEBOOKS.items():
        rebuild_script_notebook(target_name, translated_cells)
        print(f"Rebuilt {target_name} from its MATLAB example script")


if __name__ == "__main__":
    main()
