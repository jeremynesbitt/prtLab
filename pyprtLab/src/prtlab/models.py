# Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file.

"""Typed records used by the Python ray tracer."""

from __future__ import annotations

from dataclasses import dataclass, field
from html import escape
from typing import Any

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.complex128]


def _display_value(value: Any, *, limit: int | None = None) -> str:
    if isinstance(value, np.ndarray):
        text = np.array2string(
            value, precision=6, suppress_small=True, separator=", "
        )
    elif isinstance(value, dict):
        text = "{" + ", ".join(
            f"{key}: {_display_value(item)}" for key, item in value.items()
        ) + "}"
    elif isinstance(
        value, (float, complex, np.floating, np.complexfloating)
    ):
        text = f"{value:.8g}"
    else:
        text = str(value)
    if limit is not None and len(text) > limit:
        return text[: limit - 3] + "..."
    return text


def vector3() -> Array:
    return np.zeros(3, dtype=np.complex128)


def matrix3() -> Array:
    return np.eye(3, dtype=np.complex128)


@dataclass(slots=True)
class Surface:
    radius: float
    thickness: float
    surface_type: str
    surface_data: dict[str, Any] = field(default_factory=dict)
    material_type: str = "isotropic"
    index_data: dict[str, Any] = field(default_factory=dict)
    axis_data: dict[str, Any] = field(default_factory=dict)
    coating_data: dict[str, Any] = field(
        default_factory=lambda: {"type": "none"}
    )


@dataclass(slots=True)
class OpticalSystem:
    wavelength: float
    wavelength_units: str = ""
    surfaces: list[Surface] = field(default_factory=list)
    description: str = "prtLab optical system"

    def _table_rows(self) -> list[list[str]]:
        return [
            [
                str(index),
                _display_value(surface.radius),
                _display_value(surface.thickness),
                surface.surface_type,
                surface.material_type,
                _display_value(surface.index_data),
                _display_value(surface.axis_data),
                _display_value(surface.coating_data),
                _display_value(surface.surface_data),
            ]
            for index, surface in enumerate(self.surfaces, start=1)
        ]

    def to_table(self) -> str:
        """Return a compact plain-text surface table."""

        headers = [
            "#", "Radius", "Thickness", "Surface", "Material",
            "Index Data", "Axis Data", "Coating Data", "Surface Data",
        ]
        rows = [
            [_display_value(value, limit=38) for value in row]
            for row in self._table_rows()
        ]
        widths = [
            max(len(headers[column]), *(len(row[column]) for row in rows))
            if rows else len(headers[column])
            for column in range(len(headers))
        ]

        def format_row(row: list[str]) -> str:
            return " | ".join(
                value.ljust(width) for value, width in zip(row, widths, strict=True)
            )

        units = self.wavelength_units or "length units"
        lines = [
            self.description,
            f"Wavelength: {_display_value(self.wavelength)} {units}",
            format_row(headers),
            "-+-".join("-" * width for width in widths),
        ]
        lines.extend(format_row(row) for row in rows)
        return "\n".join(lines)

    def print_table(self) -> None:
        """Print the surface table for terminal or script use."""

        print(self.to_table())

    def _repr_html_(self) -> str:
        """Render the optical prescription as a table in Jupyter."""

        headers = [
            "#", "Radius", "Thickness", "Surface", "Material",
            "Index Data", "Axis Data", "Coating Data", "Surface Data",
        ]
        header_html = "".join(f"<th>{escape(item)}</th>" for item in headers)
        row_html = "".join(
            "<tr>" + "".join(
                f"<td><code>{escape(value)}</code></td>" for value in row
            ) + "</tr>"
            for row in self._table_rows()
        )
        units = self.wavelength_units or "length units"
        return (
            '<div class="prtlab-optical-system">'
            f"<strong>{escape(self.description)}</strong><br>"
            f"Wavelength: {_display_value(self.wavelength)} {escape(units)}"
            '<table style="border-collapse:collapse;margin-top:0.5em">'
            f"<thead><tr>{header_html}</tr></thead>"
            f"<tbody>{row_html}</tbody></table></div>"
        )


@dataclass(slots=True)
class TraceOptions:
    max_branches: int = 256
    min_flux: float = 100 * np.finfo(float).eps
    min_relative_flux: float = 0.0
    min_amplitude: float = 0.0
    keep_diagnostics: bool = True
    encode_propagation_phase_in_p: bool = False
    max_interactions: int = 24
    face_tolerance: float = 1e-9


@dataclass(slots=True)
class RayBranch:
    id: int = 0
    parent_id: int = 0
    surface_index: int = 0
    mode: str = ""
    branch_type: str = ""
    medium_type: str = ""
    position: Array = field(default_factory=vector3)
    k: Array = field(default_factory=vector3)
    s_direction: Array = field(default_factory=vector3)
    mode_e: Array = field(default_factory=vector3)
    mode_h: Array = field(default_factory=vector3)
    field_e: Array = field(default_factory=vector3)
    field_h: Array = field(default_factory=vector3)
    p_matrix: Array = field(default_factory=matrix3)
    q_matrix: Array = field(default_factory=matrix3)
    o_matrix: Array = field(default_factory=matrix3)
    local_basis: dict[str, Any] = field(default_factory=dict)
    amplitude: float = 0.0
    flux: float = 0.0
    opl: float = 0.0
    active: bool = True
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def S(self) -> Array:
        """MATLAB-compatible name for the energy direction."""

        return self.s_direction

    @property
    def P(self) -> Array:
        return self.p_matrix

    @property
    def Q(self) -> Array:
        return self.q_matrix

    @property
    def O(self) -> Array:
        return self.o_matrix


@dataclass(slots=True)
class Interaction:
    surface_index: int
    case_name: str
    position: Array
    normal: Array
    incident: RayBranch
    children: list[RayBranch] = field(default_factory=list)
    frames: dict[str, Any] = field(default_factory=dict)
    p_matrices: dict[str, Array] = field(default_factory=dict)
    q_matrices: dict[str, Array] = field(default_factory=dict)
    coefficients: dict[str, Any] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class TraceResult:
    system: Any
    options: TraceOptions
    rays: list[RayBranch] = field(default_factory=list)
    interactions: list[Interaction] = field(default_factory=list)
    final_ray_ids: list[int] = field(default_factory=list)


@dataclass(slots=True)
class Medium:
    material_type: str = "isotropic"
    index_data: dict[str, Any] = field(default_factory=lambda: {"n": 1.0})
    axis_data: dict[str, Any] = field(default_factory=dict)

    def as_surface(self, coating_data: dict[str, Any] | None = None) -> Surface:
        return Surface(
            radius=np.inf,
            thickness=0.0,
            surface_type="plane",
            material_type=self.material_type,
            index_data=dict(self.index_data),
            axis_data=dict(self.axis_data),
            coating_data={"type": "bare"} if coating_data is None else dict(coating_data),
        )


@dataclass(slots=True)
class SolidFace:
    name: str
    vertices: Array
    point: Array
    normal: Array
    coating_data: dict[str, Any] = field(default_factory=lambda: {"type": "bare"})
    surface_data: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class Solid:
    name: str
    wavelength: float
    wavelength_units: str
    outside: Medium
    material: Medium
    faces: list[SolidFace] = field(default_factory=list)
    geometry: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class Scene:
    wavelength: float
    outside: Medium = field(default_factory=Medium)
    wavelength_units: str = "um"
    name: str = "prtLab scene"
    solids: list[Solid] = field(default_factory=list)
    geometry: dict[str, Any] = field(default_factory=dict)
