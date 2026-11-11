"""Optional Meshers volume meshing and conversion to Microgen's VTK output."""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

import numpy as np
import pyvista as pv

if TYPE_CHECKING:
    from scipy.spatial.transform import Rotation

    from microgen.shape._types import BoundsType, Field


def generate_volume(
    field: Field,
    bounds: BoundsType | None,
    options: dict,
    *,
    default_cells=24,
    band=None,
    rotation: Rotation | None = None,
    center=(0.0, 0.0, 0.0),
) -> pv.UnstructuredGrid:
    """Convert native Tet4/Tet10 output, preserving placement and periodic pairs."""
    if bounds is None:
        raise ValueError("Bounds are required for Meshers volume meshing")
    options = dict(options)
    if "band" in options:
        raise ValueError(
            "The shape determines its solid region; band is not a meshing option"
        )
    try:
        meshers = importlib.import_module("meshers")
    except ModuleNotFoundError as error:
        if error.name != "meshers":
            raise
        raise ImportError(
            "backend='meshers' requires the optional Meshers package. "
            "Install a tested Meshers wheel or microgen[meshers]."
        ) from error

    options.setdefault("cells", default_cells)
    options.setdefault("threads", 1)
    options.setdefault("minimum_quality", 0.1)
    native = meshers.generate(field, bounds=bounds, band=band, **options)
    width = native.tetrahedra.shape[1]
    cell_type = {4: pv.CellType.TETRA, 10: pv.CellType.QUADRATIC_TETRA}[width]
    points = native.points.copy()
    lengths = np.asarray(bounds)[1::2] - np.asarray(bounds)[::2]
    shifts = np.diag(lengths)
    if rotation is not None:
        points = rotation.apply(points) + np.asarray(center)
        shifts = rotation.apply(shifts)
    grid = pv.UnstructuredGrid({cell_type: native.tetrahedra.copy()}, points)

    if len(native.surface):
        grid.field_data["meshers_boundary_faces"] = native.surface.copy()
        grid.field_data["meshers_boundary_tags"] = native.boundary_tags.copy()
    for axis, pairs, shift in zip("xyz", native.periodic_pairs, shifts, strict=True):
        if len(pairs):
            grid.field_data[f"periodic_pairs_{axis}"] = pairs.copy()
            grid.field_data[f"periodic_shift_{axis}"] = shift[None, :].copy()
    return grid
