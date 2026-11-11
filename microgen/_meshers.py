"""Optional Meshers volume meshing and conversion to Microgen's VTK output."""

from __future__ import annotations

import importlib
import importlib.metadata
import json
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
    """Generate Tet4/Tet10 without MMG, fallback, retries or relaxed guards.

    Generic Shape fields are already in world coordinates. TPMS and Spinodoid
    supply their local field and an explicit placement. Periodic pairs refer to
    returned point IDs; periodic shifts are expressed in world coordinates.
    """
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

    # These are corner metrics for either order. They are not curved Tet10
    # volumes or a replacement for its sampled-Jacobian/FEA checks.
    corners = native.points[native.tetrahedra[:, :4]]
    determinants = np.linalg.det(corners[:, 1:] - corners[:, :1])
    edge_sum = sum(
        np.sum((corners[:, i] - corners[:, j]) ** 2, axis=1)
        for i in range(4)
        for j in range(i + 1, 4)
    )
    grid.cell_data["MMGQuality"] = np.sqrt(432 * determinants**2 / edge_sum**3)
    grid.cell_data["CornerSignedVolume"] = determinants / 6
    grid.field_data["meshers_diagnostics"] = [json.dumps(native.diagnostics)]
    grid.field_data["meshers_version"] = [importlib.metadata.version("meshers")]
    if len(native.surface):
        grid.field_data["meshers_boundary_faces"] = native.surface.copy()
        grid.field_data["meshers_boundary_tags"] = native.boundary_tags.copy()
    for axis, pairs, shift in zip("xyz", native.periodic_pairs, shifts, strict=True):
        if len(pairs):
            grid.field_data[f"periodic_pairs_{axis}"] = pairs.copy()
            grid.field_data[f"periodic_shift_{axis}"] = shift[None, :].copy()
    return grid
