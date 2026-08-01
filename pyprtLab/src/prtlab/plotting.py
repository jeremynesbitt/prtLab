# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Matplotlib visualization for optical systems and polarized ray trees."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from .geometry import surface_sag
from .models import OpticalSystem, RayBranch, Scene, Solid, TraceResult


def polarization_ellipse_points_2d(
    field,
    *,
    center=(0.0, 0.0),
    scale=1.0,
    normalize_field=False,
    samples=101,
):
    """Return the time-domain polarization ellipse for a Jones vector."""

    value = np.asarray(field, dtype=np.complex128).reshape(-1)
    if value.size != 2:
        raise ValueError("field must be a two-element Jones vector")
    magnitude = np.linalg.norm(value)
    if magnitude <= np.finfo(float).eps:
        return np.repeat(np.asarray(center, dtype=float)[:, None], samples, axis=1)
    if normalize_field:
        value = value / magnitude
    theta = np.linspace(0, 2 * np.pi, samples)
    return (
        np.asarray(center, dtype=float).reshape(2, 1)
        + scale * np.real(np.exp(-1j * theta)[None, :] * value[:, None])
    )


def plot_jones_vector(
    field,
    *,
    center=(0.0, 0.0),
    scale=1.0,
    normalize_field=False,
    samples=101,
    color=(0.85, 0.10, 0.10),
    line_width=2.0,
    show_arrow=True,
    ax=None,
):
    """Plot a Jones vector as a time-domain polarization ellipse."""

    if ax is None:
        _, ax = plt.subplots()
    points = polarization_ellipse_points_2d(
        field, center=center, scale=scale, normalize_field=normalize_field,
        samples=samples,
    )
    artists = list(ax.plot(points[0], points[1], color=color, linewidth=line_width))
    if show_arrow and samples >= 3:
        index = max(1, round(samples * 5 / 12))
        tangent = points[:, (index + 1) % samples] - points[:, index - 1]
        if np.linalg.norm(tangent) > np.finfo(float).eps:
            tangent = tangent / np.linalg.norm(tangent)
            artists.append(ax.quiver(
                *points[:, index], *(0.18 * scale * tangent),
                color=color, angles="xy", scale_units="xy", scale=1,
            ))
    ax.set_aspect("equal", adjustable="box")
    return ax, artists


def prepare_jones_pupil(jones, mask=None, *, relative_phase_threshold=1e-4):
    """Prepare signed amplitude and sign-absorbed phase pupil arrays."""

    value = np.asarray(jones, dtype=np.complex128)
    if value.ndim != 4 or value.shape[-2:] != (2, 2):
        raise ValueError("jones must have shape (rows, columns, 2, 2)")
    rows, columns = value.shape[:2]
    if mask is None:
        x_values = np.linspace(-1, 1, columns)
        y_values = np.linspace(-1, 1, rows)
        xx, yy = np.meshgrid(x_values, y_values)
        mask = np.hypot(xx, yy) < 1
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != (rows, columns):
        raise ValueError("mask must match the first two Jones-pupil dimensions")

    signed_amplitude = np.full(value.shape, np.nan, dtype=float)
    phase = np.full(value.shape, np.nan, dtype=float)
    phase_mask = np.zeros(value.shape, dtype=bool)
    phase_offset = np.full((2, 2), np.nan)
    for row in range(2):
        for column in range(2):
            element = value[:, :, row, column]
            available = np.abs(element[mask])
            if available.size == 0:
                continue
            threshold = relative_phase_threshold * np.max(available)
            valid = mask & (np.abs(element) > threshold)
            if np.any(valid):
                offset = np.angle(np.mean(element[valid]))
            else:
                offset = 0.0
            rotated = element * np.exp(-1j * offset)
            sign = np.sign(np.real(rotated))
            sign[sign == 0] = 1
            amplitude = np.real(rotated)
            amplitude[~mask] = np.nan
            element_phase = np.angle(rotated * sign)
            element_phase[~valid] = np.nan
            signed_amplitude[:, :, row, column] = amplitude
            phase[:, :, row, column] = element_phase
            phase_mask[:, :, row, column] = valid
            phase_offset[row, column] = offset
    return {
        "signed_amplitude": signed_amplitude,
        "phase": phase,
        "phase_mask": phase_mask,
        "phase_offset": phase_offset,
        "pupil_mask": mask,
    }


def plot_jones_pupil(jones, mask=None, *, relative_phase_threshold=1e-4):
    """Plot signed-amplitude and sign-absorbed phase Jones pupils."""

    data = prepare_jones_pupil(
        jones, mask, relative_phase_threshold=relative_phase_threshold
    )
    amplitude_figure, amplitude_axes = plt.subplots(2, 2, squeeze=False)
    phase_figure, phase_axes = plt.subplots(2, 2, squeeze=False)
    for row in range(2):
        for column in range(2):
            amplitude_image = amplitude_axes[row, column].imshow(
                np.flipud(data["signed_amplitude"][:, :, row, column])
            )
            amplitude_figure.colorbar(amplitude_image, ax=amplitude_axes[row, column])
            phase_image = phase_axes[row, column].imshow(
                np.flipud(data["phase"][:, :, row, column])
            )
            phase_figure.colorbar(phase_image, ax=phase_axes[row, column])
            amplitude_axes[row, column].set_title(f"J{row + 1}{column + 1}")
            phase_axes[row, column].set_title(f"J{row + 1}{column + 1}")
    amplitude_figure.suptitle("Signed Amplitude of Jones Pupil")
    phase_figure.suptitle("Phase of Jones Pupil [rad]")
    return (amplitude_figure, amplitude_axes), (phase_figure, phase_axes), data


def plot_polarization_ellipses_across_field(
    x, y, jones, input_field, mask=None, *, scale=None, ax=None
):
    """Plot output polarization ellipses across a sampled field or pupil."""

    xx, yy = np.asarray(x), np.asarray(y)
    if xx.shape != yy.shape:
        raise ValueError("x and y must have matching shapes")
    matrix = np.asarray(jones, dtype=np.complex128)
    if matrix.shape == (*xx.shape, 4):
        matrix = matrix.reshape(*xx.shape, 2, 2)
    if matrix.shape != (*xx.shape, 2, 2):
        raise ValueError("jones must have shape x.shape + (2, 2) or x.shape + (4,)")
    if mask is None:
        mask = np.ones(xx.shape, dtype=bool)
    mask = np.asarray(mask, dtype=bool)
    if scale is None:
        unique_y = np.unique(yy)
        scale = 0.5 * np.min(np.diff(unique_y)) if unique_y.size > 1 else 0.5
    if ax is None:
        _, ax = plt.subplots()
    artists = []
    input_value = np.asarray(input_field, dtype=np.complex128).reshape(2)
    for index in np.ndindex(xx.shape):
        if not mask[index]:
            continue
        output = matrix[index] @ input_value
        _, handles = plot_jones_vector(
            output, center=(-xx[index], yy[index]), scale=scale,
            normalize_field=True, show_arrow=False, ax=ax,
        )
        artists.extend(handles)
    ax.set_title("Output Polarization")
    return ax, artists


def plot_basis_vectors_on_sphere(
    k_grid,
    x_local,
    y_local=None,
    *,
    show_sphere=True,
    x_color="r",
    y_color=(0, 0.7, 0),
    arrow_scale=0.15,
    sphere_alpha=0.15,
    sphere_color=(0.7, 0.7, 1),
    sphere_samples=30,
    view=(-37.5, 30),
    line_width=1.0,
    ax=None,
):
    """Plot local transverse basis vectors on the unit propagation sphere."""

    k_values = np.asarray(k_grid, dtype=float).reshape(-1, 3)
    if k_values.shape[0] == 0:
        raise ValueError("k_grid must contain at least one propagation vector")
    x_values = (
        None if x_local is None else np.asarray(x_local, dtype=float).reshape(-1, 3)
    )
    y_values = (
        None if y_local is None else np.asarray(y_local, dtype=float).reshape(-1, 3)
    )
    if x_values is not None and x_values.shape != k_values.shape:
        raise ValueError("x_local must have the same shape as k_grid")
    if y_values is not None and y_values.shape != k_values.shape:
        raise ValueError("y_local must have the same shape as k_grid")
    if ax is None:
        _, ax = plt.subplots(subplot_kw={"projection": "3d"})
    handles = {"sphere": None, "x": None, "y": None}
    if show_sphere:
        azimuth = np.linspace(0, 2 * np.pi, sphere_samples + 1)
        elevation = np.linspace(-np.pi / 2, np.pi / 2, sphere_samples + 1)
        aa, ee = np.meshgrid(azimuth, elevation)
        handles["sphere"] = ax.plot_surface(
            np.cos(ee) * np.cos(aa),
            np.cos(ee) * np.sin(aa),
            np.sin(ee),
            color=sphere_color,
            alpha=sphere_alpha,
            linewidth=0.15,
        )
    if x_values is not None:
        handles["x"] = ax.quiver(
            *k_values.T,
            *(arrow_scale * x_values).T,
            color=x_color,
            linewidth=line_width,
            normalize=False,
        )
    if y_values is not None:
        handles["y"] = ax.quiver(
            *k_values.T,
            *(arrow_scale * y_values).T,
            color=y_color,
            linewidth=line_width,
            normalize=False,
        )
    ax.set(xlabel="x", ylabel="y", zlabel="z")
    ax.set_box_aspect((1, 1, 1))
    view_values = np.asarray(view, dtype=float).reshape(-1)
    if view_values.size == 2:
        azimuth, elevation = view_values
    elif view_values.size == 3:
        if np.linalg.norm(view_values) == 0:
            raise ValueError("a 3D view direction must be nonzero")
        azimuth = np.rad2deg(np.arctan2(view_values[1], view_values[0]))
        elevation = np.rad2deg(
            np.arctan2(view_values[2], np.hypot(*view_values[:2]))
        )
    else:
        raise ValueError("view must contain two angles or a 3D direction")
    ax.view_init(elev=elevation, azim=azimuth)
    ax.grid(True)
    return ax, handles


def _vertex_positions(system: OpticalSystem) -> np.ndarray:
    return np.cumsum([surface.thickness for surface in system.surfaces])


def _is_lens_pair(system: OpticalSystem, index: int) -> bool:
    next_surface = system.surfaces[index + 1]
    if next_surface.material_type != "isotropic":
        return False
    index_value = complex(next_surface.index_data["n"])
    return abs(index_value.imag) <= np.finfo(float).eps and index_value.real != 1


def plot_prt_lens_cross_section(
    system: OpticalSystem,
    *,
    num_points=120,
    lens_color=(0.60, 0.80, 1.00),
    filter_color=(0.60, 1.00, 0.75),
    default_clear_aperture=1.0,
    show_labels=True,
    show_stop=False,
    stop_z=0.20,
    entrance_pupil_radius=5.57 / (2 * 2.8),
    ax=None,
):
    """Plot a sequential optical system in the y-z meridional plane."""

    if ax is None:
        _, ax = plt.subplots()
    count = max(0, len(system.surfaces) - 1)
    vertices = _vertex_positions(system)
    handles = {"elements": [], "surfaces": [], "labels": [], "stop": []}
    index = 0
    while index < count:
        surface = system.surfaces[index]
        paired = index < count - 1
        if paired and _is_lens_pair(system, index):
            second = system.surfaces[index + 1]
            aperture = min(
                _clear_aperture(surface, default_clear_aperture),
                _clear_aperture(second, default_clear_aperture),
            )
            y_values = np.linspace(-aperture, aperture, 2 * num_points - 1)
            z_first = np.array(
                [vertices[index] + surface_sag(0, y, surface) for y in y_values]
            )
            z_second = np.array(
                [
                    vertices[index + 1] + surface_sag(0, y, second)
                    for y in y_values
                ]
            )
            handles["elements"].append(
                ax.fill(
                    np.r_[z_first, z_second[::-1]],
                    np.r_[y_values, y_values[::-1]],
                    color=lens_color,
                    edgecolor="k",
                    linewidth=1.2,
                )[0]
            )
            index += 2
            continue
        if (
            paired
            and surface.surface_type == "plane"
            and system.surfaces[index + 1].surface_type == "plane"
        ):
            second = system.surfaces[index + 1]
            aperture = min(
                _clear_aperture(surface, default_clear_aperture),
                _clear_aperture(second, default_clear_aperture),
            )
            handles["elements"].append(
                ax.fill(
                    [vertices[index], vertices[index + 1],
                     vertices[index + 1], vertices[index]],
                    [-aperture, -aperture, aperture, aperture],
                    color=filter_color,
                    edgecolor="k",
                    linewidth=1.2,
                )[0]
            )
            index += 2
            continue
        aperture = _clear_aperture(surface, default_clear_aperture)
        y_values = np.linspace(-aperture, aperture, num_points)
        z_values = np.array(
            [vertices[index] + surface_sag(0, y, surface) for y in y_values]
        )
        handles["surfaces"].extend(
            ax.plot(z_values, y_values, color="k", linewidth=1.2)
        )
        index += 1

    if show_stop:
        extension = 0.25
        handles["stop"].extend(
            ax.plot(
                [stop_z, stop_z],
                [entrance_pupil_radius, entrance_pupil_radius + extension],
                color="k",
                linewidth=3,
            )
        )
        handles["stop"].extend(
            ax.plot(
                [stop_z, stop_z],
                [-entrance_pupil_radius - extension, -entrance_pupil_radius],
                color="k",
                linewidth=3,
            )
        )
    if count:
        z_min = float(np.min(vertices[:count]) - 0.5)
        z_max = float(np.max(vertices[:count]) + 0.7)
    else:
        z_min, z_max = -0.5, 0.7
    handles["axis"] = ax.plot(
        [z_min, z_max], [0, 0], color="k", linestyle=":", linewidth=0.8
    )
    if show_labels:
        for surface_index in range(count):
            surface = system.surfaces[surface_index]
            aperture = _clear_aperture(surface, default_clear_aperture)
            z_label = vertices[surface_index] + surface_sag(
                0, 0.85 * aperture, surface
            )
            handles["labels"].append(
                ax.text(
                    z_label,
                    aperture + 0.12,
                    str(surface_index + 1),
                    fontsize=9,
                    horizontalalignment="center",
                    fontweight="bold",
                )
            )
    units = system.wavelength_units or "length units"
    ax.set(
        xlabel=f"z ({units})",
        ylabel=f"y ({units})",
        xlim=(z_min, z_max),
    )
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    return ax, handles


def plot_prt_ray_trace(
    ray_outputs,
    *,
    color=(0.85, 0.10, 0.10),
    line_width=1.5,
    pre_extend=0.3,
    post_extend=0.5,
    plot_incomplete=True,
    incomplete_line_style="--",
    auto_expand_axes=True,
    axis_padding_fraction=0.04,
    color_by_mode=True,
    mode_colors=None,
    ax=None,
):
    """Overlay sequential ray histories in the y-z meridional plane."""

    if ax is None:
        _, ax = plt.subplots()
    outputs = (
        list(ray_outputs)
        if isinstance(ray_outputs, (list, tuple))
        else [ray_outputs]
    )
    default_colors = PlotOptions().mode_colors
    colors_by_mode = default_colors if mode_colors is None else mode_colors
    handles = []
    plotted = []
    for result in outputs:
        final_ids = result.final_ray_ids
        if not final_ids and plot_incomplete and result.rays:
            positions = [np.real(result.rays[0].position)]
            positions.extend(np.real(item.position) for item in result.interactions)
            points = np.asarray(positions)
            pre = points[0] - pre_extend * np.real(result.rays[0].S)
            handles.extend(
                ax.plot(
                    [pre[2], points[0, 2]],
                    [pre[1], points[0, 1]],
                    linestyle=incomplete_line_style,
                    color=color,
                    linewidth=line_width,
                )
            )
            handles.extend(
                ax.plot(
                    points[:, 2],
                    points[:, 1],
                    linestyle=incomplete_line_style,
                    color=color,
                    linewidth=line_width,
                )
            )
            plotted.extend([pre, *points])
            continue
        for final_id in final_ids:
            history = _history(result, final_id)
            rays = [result.rays[ray_id - 1] for ray_id in history]
            positions = np.asarray([np.real(ray.position) for ray in rays])
            pre = positions[0] - pre_extend * np.real(rays[0].S)
            initial_color = (
                _ray_color(rays[0], PlotOptions(mode_colors=colors_by_mode))
                if color_by_mode
                else color
            )
            handles.extend(
                ax.plot(
                    [pre[2], positions[0, 2]],
                    [pre[1], positions[0, 1]],
                    color=initial_color,
                    linewidth=line_width,
                )
            )
            for segment_index in range(1, len(rays)):
                ray = rays[segment_index - 1]
                segment_color = (
                    _ray_color(ray, PlotOptions(mode_colors=colors_by_mode))
                    if color_by_mode
                    else color
                )
                handles.extend(
                    ax.plot(
                        positions[segment_index - 1:segment_index + 1, 2],
                        positions[segment_index - 1:segment_index + 1, 1],
                        color=segment_color,
                        linewidth=line_width,
                    )
                )
            final_ray = rays[-1]
            post = positions[-1] + post_extend * np.real(final_ray.S)
            final_color = (
                _ray_color(final_ray, PlotOptions(mode_colors=colors_by_mode))
                if color_by_mode
                else color
            )
            handles.extend(
                ax.plot(
                    [positions[-1, 2], post[2]],
                    [positions[-1, 1], post[1]],
                    color=final_color,
                    linewidth=line_width,
                )
            )
            plotted.extend([pre, *positions, post])
    if auto_expand_axes and plotted:
        points = np.asarray(plotted)[:, [2, 1]]
        x_limits, y_limits = ax.get_xlim(), ax.get_ylim()
        low = np.minimum([x_limits[0], y_limits[0]], np.min(points, axis=0))
        high = np.maximum([x_limits[1], y_limits[1]], np.max(points, axis=0))
        padding = axis_padding_fraction * np.maximum(
            high - low, np.finfo(float).eps
        )
        ax.set_xlim(low[0] - padding[0], high[0] + padding[0])
        ax.set_ylim(low[1] - padding[1], high[1] + padding[1])
    return ax, handles


@dataclass(slots=True)
class PlotOptions:
    show_system: bool = True
    show_rays: bool = True
    show_polarization: bool = True
    ray_selection: str = "all"
    polarization_at: str = "interfaces"
    field_name: str = "field_e"
    polarization_scale: float | str = "auto"
    polarization_scale_fraction: float = 0.10
    normalize_polarization: bool = True
    samples_per_ellipse: int = 101
    surface_samples: int = 36
    default_clear_aperture: float = 1.0
    surface_alpha: float = 0.18
    line_width: float = 1.5
    glyph_line_width: float = 1.4
    pre_extend: float = 0.0
    post_extend: float = 0.0
    axis_padding_fraction: float = 0.16
    equal_axes: bool = True
    show_legend: bool = False
    color_by_mode: bool = True
    mode_colors: dict[str, tuple[float, float, float]] = field(
        default_factory=lambda: {
            "input": (0.12, 0.12, 0.12),
            "isotropic": (0.12, 0.12, 0.12),
            "transmitted": (0.12, 0.12, 0.12),
            "ordinary": (0.00, 0.30, 0.90),
            "extraordinary": (0.00, 0.55, 0.20),
            "reflected": (0.95, 0.45, 0.05),
        }
    )
    polarization_color: tuple[float, float, float] = (0.85, 0.10, 0.10)
    view: tuple[float, float] = (20, 35)


def _options(value: PlotOptions | dict[str, Any] | None) -> PlotOptions:
    if value is None:
        return PlotOptions()
    if isinstance(value, PlotOptions):
        return value
    aliases = {
        "ShowSystem": "show_system", "ShowRays": "show_rays",
        "ShowPolarization": "show_polarization", "RaySelection": "ray_selection",
        "PolarizationAt": "polarization_at", "FieldName": "field_name",
        "PolarizationScale": "polarization_scale",
        "PolarizationScaleFraction": "polarization_scale_fraction",
        "NormalizePolarization": "normalize_polarization",
        "SamplesPerEllipse": "samples_per_ellipse", "SurfaceSamples": "surface_samples",
        "DefaultClearAperture": "default_clear_aperture", "SurfaceAlpha": "surface_alpha",
        "LineWidth": "line_width", "GlyphLineWidth": "glyph_line_width",
        "PreExtend": "pre_extend", "PostExtend": "post_extend",
        "AxisPaddingFraction": "axis_padding_fraction", "EqualAxes": "equal_axes",
        "ShowLegend": "show_legend", "ColorByMode": "color_by_mode", "View": "view",
    }
    return PlotOptions(**{aliases.get(key, key): item for key, item in value.items()})


def _ray_color(ray: RayBranch, options: PlotOptions):
    if not options.color_by_mode:
        return options.mode_colors["isotropic"]
    key = ray.mode
    if ray.branch_type in {"input", "reflected"} or key not in options.mode_colors:
        key = ray.branch_type
    return options.mode_colors.get(key, options.mode_colors["isotropic"])


def _history(result: TraceResult, final_id: int) -> list[int]:
    ray = result.rays[final_id - 1]
    if "history" in ray.metadata:
        return list(ray.metadata["history"])
    history = []
    while ray.id:
        history.append(ray.id)
        if ray.parent_id == 0:
            break
        ray = result.rays[ray.parent_id - 1]
    return history[::-1]


def _selected_final_ids(result: TraceResult, selection: str) -> list[int]:
    if selection == "all":
        return result.final_ray_ids
    if selection != "dominant":
        raise ValueError("ray_selection must be 'all' or 'dominant'")
    return [max(result.final_ray_ids, key=lambda ray_id: abs(result.rays[ray_id - 1].flux))] if result.final_ray_ids else []


def _plot_solid(ax, solid: Solid, options: PlotOptions):
    color = (0.60, 0.80, 1.00) if solid.material.material_type == "isotropic" else (0.85, 0.85, 0.85)
    polygons = [np.real(face.vertices.T) for face in solid.faces]
    collection = Poly3DCollection(
        polygons, facecolor=color, edgecolor=(0, 0, 0, 0.35),
        linewidth=0.8, alpha=options.surface_alpha,
    )
    ax.add_collection3d(collection)
    return np.vstack(polygons), [collection]


def _clear_aperture(surface, default: float) -> float:
    for key in ("clearAperture", "ClearAperture", "semiDiameter"):
        if key in surface.surface_data:
            return float(surface.surface_data[key])
    return default


def _plot_optical_system(ax, system: OpticalSystem, options: PlotOptions):
    points = []
    artists = []
    vertex_z = 0.0
    phi = np.linspace(0, 2 * np.pi, options.surface_samples)
    radial = np.linspace(0, 1, options.surface_samples)
    rr, pp = np.meshgrid(radial, phi)
    for surface in system.surfaces[:-1]:
        vertex_z += surface.thickness
        aperture = _clear_aperture(surface, options.default_clear_aperture)
        xx = aperture * rr * np.cos(pp)
        yy = aperture * rr * np.sin(pp)
        zz = np.empty_like(xx)
        for index in np.ndindex(xx.shape):
            zz[index] = vertex_z + surface_sag(xx[index], yy[index], surface)
        artist = ax.plot_surface(
            xx, yy, zz, color=(0.60, 0.80, 1.00),
            alpha=options.surface_alpha, linewidth=0.1,
        )
        artists.append(artist)
        points.append(np.column_stack((xx.ravel(), yy.ravel(), zz.ravel())))
    return (np.vstack(points) if points else np.empty((0, 3))), artists


def _plot_geometry(ax, system, options: PlotOptions):
    if isinstance(system, Scene):
        chunks, artists = [], []
        for solid in system.solids:
            points, handles = _plot_solid(ax, solid, options)
            chunks.append(points)
            artists.extend(handles)
        return (np.vstack(chunks) if chunks else np.empty((0, 3))), artists
    if isinstance(system, Solid):
        return _plot_solid(ax, system, options)
    if isinstance(system, OpticalSystem):
        return _plot_optical_system(ax, system, options)
    raise TypeError("system must be an OpticalSystem, Solid, or Scene")


def _plot_ray_history(ax, result: TraceResult, history: list[int], options: PlotOptions):
    rays = [result.rays[ray_id - 1] for ray_id in history]
    positions = np.array([np.real(ray.position) for ray in rays])
    points = [positions]
    artists = []
    if options.pre_extend > 0:
        pre = positions[0] - options.pre_extend * np.real(rays[0].S)
        artists.extend(ax.plot(*np.column_stack((pre, positions[0])), color=_ray_color(rays[0], options), linewidth=options.line_width))
        points.append(pre[None, :])
    for index in range(1, len(rays)):
        segment = positions[index - 1:index + 1].T
        artists.extend(ax.plot(*segment, color=_ray_color(rays[index - 1], options), linewidth=options.line_width))
    if options.post_extend > 0:
        post = positions[-1] + options.post_extend * np.real(rays[-1].S)
        artists.extend(ax.plot(*np.column_stack((positions[-1], post)), color=_ray_color(rays[-1], options), linewidth=options.line_width))
        points.append(post[None, :])
    return np.vstack(points), artists


def plot_polarization_glyph_3d(ax, center, field, *, scale=1.0, color=(0.85, 0.10, 0.10), normalize_field=True, samples=101, line_width=1.4):
    """Plot the real time-domain locus of a complex 3D electric field."""
    field = np.asarray(field, dtype=np.complex128).reshape(3)
    magnitude = np.linalg.norm(field)
    if magnitude <= np.finfo(float).eps:
        return []
    if normalize_field:
        field = field / magnitude
    field *= scale
    center = np.asarray(center, dtype=float).reshape(3)
    theta = np.linspace(0, 2 * np.pi, samples)
    points = center[:, None] + np.real(np.exp(-1j * theta)[None, :] * field[:, None])
    artists = list(ax.plot(*points, color=color, linewidth=line_width))
    theta0, delta = -np.pi / 6, 1e-2
    p0 = center + np.real(np.exp(-1j * theta0) * field)
    p1 = center + np.real(np.exp(-1j * (theta0 + delta)) * field)
    tangent = p1 - p0
    if np.linalg.norm(tangent) > np.finfo(float).eps:
        tangent /= np.linalg.norm(tangent)
        artists.append(ax.quiver(*p0, *(0.35 * scale * tangent), color=color, linewidth=1.0))
    return artists


def _set_equal_axes(ax, points: np.ndarray, padding_fraction: float):
    low, high = np.min(points, axis=0), np.max(points, axis=0)
    center = (low + high) / 2
    radius = max(np.max(high - low) / 2, np.finfo(float).eps)
    radius *= 1 + padding_fraction
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)


def plot_prt_system_3d(system, ray_outputs, options: PlotOptions | dict[str, Any] | None = None, *, ax=None):
    """Plot system geometry, branched ray paths, and polarization ellipses."""
    options = _options(options)
    if ax is None:
        _, ax = plt.subplots(subplot_kw={"projection": "3d"})
    if ray_outputs is None:
        outputs = []
    elif isinstance(ray_outputs, (list, tuple)):
        outputs = list(ray_outputs)
    else:
        outputs = [ray_outputs]
    all_points = []
    handles = {"surfaces": [], "rays": [], "glyphs": []}
    if options.show_system:
        points, artists = _plot_geometry(ax, system, options)
        all_points.append(points)
        handles["surfaces"].extend(artists)

    extent_points = [chunk for chunk in all_points if chunk.size]
    for result in outputs:
        for ray_id in _selected_final_ids(result, options.ray_selection):
            history = _history(result, ray_id)
            positions = np.array([np.real(result.rays[index - 1].position) for index in history])
            extent_points.append(positions)
    if isinstance(options.polarization_scale, str):
        if options.polarization_scale != "auto":
            raise ValueError("polarization_scale must be 'auto' or positive")
        combined = np.vstack(extent_points) if extent_points else np.zeros((1, 3))
        scale = max(options.polarization_scale_fraction * np.max(np.ptp(combined, axis=0)), np.finfo(float).eps)
    else:
        scale = float(options.polarization_scale)
        if scale <= 0:
            raise ValueError("polarization_scale must be positive")

    plotted_glyphs = set()
    for result in outputs:
        for ray_id in _selected_final_ids(result, options.ray_selection):
            history = _history(result, ray_id)
            if options.show_rays:
                points, artists = _plot_ray_history(ax, result, history, options)
                all_points.append(points)
                handles["rays"].extend(artists)
            if options.show_polarization:
                glyph_ids = [history[-1]] if options.polarization_at == "final" else history
                for glyph_id in glyph_ids:
                    ray = result.rays[glyph_id - 1]
                    key = (*np.round(np.real(ray.position), 12), ray.mode)
                    if options.polarization_at == "interfaces" and key in plotted_glyphs:
                        continue
                    plotted_glyphs.add(key)
                    field = getattr(ray, options.field_name)
                    handles["glyphs"].extend(plot_polarization_glyph_3d(
                        ax, np.real(ray.position), field, scale=scale,
                        color=options.polarization_color,
                        normalize_field=options.normalize_polarization,
                        samples=options.samples_per_ellipse,
                        line_width=options.glyph_line_width,
                    ))

    units = getattr(system, "wavelength_units", "") or "length units"
    ax.set(xlabel=f"x ({units})", ylabel=f"y ({units})", zlabel=f"z ({units})")
    ax.view_init(elev=options.view[0], azim=options.view[1])
    ax.grid(True)
    combined = np.vstack([chunk for chunk in all_points if chunk.size]) if any(chunk.size for chunk in all_points) else np.empty((0, 3))
    if combined.size:
        if options.equal_axes:
            _set_equal_axes(ax, combined, options.axis_padding_fraction)
        else:
            low, high = np.min(combined, axis=0), np.max(combined, axis=0)
            padding = options.axis_padding_fraction * np.maximum(high - low, np.finfo(float).eps)
            ax.set_xlim(low[0] - padding[0], high[0] + padding[0])
            ax.set_ylim(low[1] - padding[1], high[1] + padding[1])
            ax.set_zlim(low[2] - padding[2], high[2] + padding[2])
    if options.show_legend:
        entries = [Line2D([0], [0], color=color, label=name, linewidth=options.line_width) for name, color in options.mode_colors.items() if name in {"input", "isotropic", "ordinary", "extraordinary", "reflected"}]
        entries.append(Line2D([0], [0], color=options.polarization_color, label="polarization", linewidth=options.glyph_line_width))
        ax.legend(handles=entries, loc="best")
    return ax, handles
