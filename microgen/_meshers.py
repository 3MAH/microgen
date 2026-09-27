"""Mesh bounded negative implicit fields with the meshers 0.1.0 backend."""

from __future__ import annotations

import meshers
import numpy as np
import pyvista as pv


def generate(
    field, bounds, resolution=50, *, periodic=(False, False, False), **options
):
    """Use grid-point resolution, with half the largest cell spacing as tolerance.

    ``meshers`` options, including a stricter ``geometry_tolerance``, are passed
    through. Callbacks are explicit because microgen fields can call VTK,
    interpolation, and automatic differentiation code outside the compiler.
    """
    counts = np.asarray(resolution)
    if counts.ndim == 0:
        counts = np.repeat(counts, 3)
    if (
        counts.shape != (3,)
        or not np.issubdtype(counts.dtype, np.integer)
        or np.any(counts < 3)
    ):
        raise ValueError("resolution must contain at least three grid points per axis")
    cells = counts - 1
    bounds = np.asarray(bounds, dtype=float)
    if (
        bounds.shape != (6,)
        or not np.all(np.isfinite(bounds))
        or np.any(bounds[1::2] <= bounds[::2])
    ):
        raise ValueError("bounds must contain three finite increasing intervals")
    options.setdefault(
        "geometry_tolerance", float(np.max((bounds[1::2] - bounds[::2]) / cells) / 2)
    )
    options.setdefault("compile", False)

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
    for axis, pairs in zip("xyz", result.periodic_pairs, strict=True):
        if len(pairs):
            grid.field_data[f"periodic_pairs_{axis}"] = pairs.copy()
    if result.periodic_transforms is not None:
        grid.field_data["periodic_transforms"] = result.periodic_transforms.reshape(
            3, 16
        ).copy()
    return grid


def surface(result):
    """Return the closed solid boundary with face tags and original node IDs."""
    if not len(result.surface):
        return pv.PolyData()
    faces = np.column_stack((np.full(len(result.surface), 3), result.surface))
    mesh = pv.PolyData(result.points.copy(), faces)
    mesh.cell_data["BoundaryTag"] = result.boundary_tags.copy()
    mesh.point_data["meshers_node_id"] = np.arange(len(result.points))
    # Remove unused interior points. Do not merge distinct boundary vertices.
    mesh = mesh.clean(point_merging=False)
    ids = np.full(len(result.points), -1, dtype=np.int64)
    ids[mesh.point_data["meshers_node_id"]] = np.arange(mesh.n_points)
    for axis, pairs in zip("xyz", result.periodic_pairs, strict=True):
        if len(pairs):
            mesh.field_data[f"periodic_pairs_{axis}"] = ids[pairs]
    if result.periodic_transforms is not None:
        mesh.field_data["periodic_transforms"] = result.periodic_transforms.reshape(
            3, 16
        ).copy()
    return mesh


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
