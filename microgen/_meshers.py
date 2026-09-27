"""Mesh bounded negative implicit fields with the meshers 0.1.0 backend."""

from __future__ import annotations

import meshers
import json
import numpy as np
import pyvista as pv


def interpolate_offset(points, values):
    """Interpolate fixed nodal offsets without batch-dependent normalization."""
    from scipy.interpolate import (
        LinearNDInterpolator,
        NearestNDInterpolator,
        RegularGridInterpolator,
    )

    axes = [np.unique(points[:, axis]) for axis in range(3)]
    if np.prod([len(axis) for axis in axes]) == len(points):
        data = np.empty(tuple(len(axis) for axis in axes))
        indices = tuple(
            np.searchsorted(axis, points[:, i]) for i, axis in enumerate(axes)
        )
        data[indices] = values
        interpolation = RegularGridInterpolator(
            axes, data, bounds_error=False, fill_value=None
        )

        def field(x, y, z):
            return interpolation(np.column_stack((x, y, z)))
    else:
        linear = LinearNDInterpolator(points, values)
        nearest = NearestNDInterpolator(points, values)

        def field(x, y, z):
            result = linear(x, y, z)
            missing = np.isnan(result)
            result[missing] = nearest(x[missing], y[missing], z[missing])
            return result

    return field


def generate(
    field, bounds, resolution=50, *, periodic=(False, False, False), **options
):
    """Convert microgen grid-point resolution to meshers background cells.

    ``meshers`` options, including a stricter ``geometry_tolerance``, are passed
    through. Supported fields compile inside the meshers wheel; VTK and
    interpolated fields use the callback path.
    """
    counts = np.asarray(resolution)
    if counts.ndim == 0:
        counts = np.repeat(counts, 3)
    if (
        counts.shape != (3,)
        or not np.issubdtype(counts.dtype, np.integer)
        or np.any(counts < 5)
        or np.any(counts > 129)
    ):
        raise ValueError(
            "meshers requires 5 to 129 grid points per axis, including repeats"
        )
    cells = counts - 1
    bounds = np.asarray(bounds, dtype=float)
    if (
        bounds.shape != (6,)
        or not np.all(np.isfinite(bounds))
        or np.any(bounds[1::2] <= bounds[::2])
    ):
        raise ValueError("bounds must contain three finite increasing intervals")
    options.setdefault("geometry_tolerance", 0.01)
    options.setdefault("compile", True)

    def callback(x, y, z):
        return np.broadcast_to(
            np.asarray(field(x, y, z), dtype=np.float64), x.shape
        ).copy()

    if isinstance(field, dict):
        return meshers.generate_intersection(
            field, bounds=bounds, cells=cells, periodic=periodic, **options
        )
    return meshers.generate(
        field if options["compile"] else callback,
        bounds=bounds,
        cells=cells,
        periodic=periodic,
        **options,
    )


def volume(result):
    """Copy native arrays into a mutable PyVista tetrahedral mesh."""
    grid = pv.UnstructuredGrid(
        {pv.CellType.TETRA: result.tetrahedra.copy()}, result.points.copy()
    )
    p = result.points[result.tetrahedra]
    determinants = np.linalg.det(p[:, 1:] - p[:, :1])
    edge_sum = sum(
        np.sum((p[:, i] - p[:, j]) ** 2, axis=1)
        for i in range(4)
        for j in range(i + 1, 4)
    )
    grid.cell_data["Volume"] = determinants / 6
    grid.cell_data["MMGQuality"] = np.sqrt(432 * determinants**2 / edge_sum**3)
    grid.field_data["meshers_diagnostics"] = [json.dumps(result.diagnostics)]
    for axis, pairs in zip("xyz", result.periodic_pairs, strict=True):
        if len(pairs):
            grid.field_data[f"periodic_pairs_{axis}"] = pairs.copy()
    if result.periodic_transforms is not None:
        grid.field_data["periodic_transforms"] = result.periodic_transforms.reshape(
            3, 16
        ).copy()
    return grid


def transform(result, rotation, center):
    """Apply a rigid placement to points and any periodic transformations."""
    result.points[:] = rotation.apply(result.points) + np.asarray(center)
    if result.periodic_transforms is not None:
        placement = np.eye(4)
        placement[:3, :3] = rotation.as_matrix()
        placement[:3, 3] = center
        result.periodic_transforms[:] = (
            placement @ result.periodic_transforms @ np.linalg.inv(placement)
        )
    return result
