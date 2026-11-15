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
    constraints=None,
    coordinate_map=None,
    periodic_transforms=None,
    weld_axis=None,
    adaptive_field=None,
    inverse_map=None,
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

    if (
        constraints is not None
        and options.get("adapt", False)
        and adaptive_field is not None
    ):
        seeded = (
            getattr(meshers, "ADAPTIVE_MAPPED_INTERSECTIONS", False)
            and inverse_map is not None
            and options.get("order", 1) == 1
            and options.get("threads", 1) == 1
            and not options.get("curved_fit", False)
            and not options.get("interior_refine_nodes", 0)
            and not any(key in options for key in ("gradient", "size_field", "layers"))
        )
        if not seeded:
            field, field_bounds = adaptive_field
            return generate_volume(
                field,
                field_bounds,
                options,
                default_cells=default_cells,
                rotation=rotation,
                center=center,
            )
    options.setdefault("cells", default_cells)
    options.setdefault("threads", 1)
    options.setdefault("minimum_quality", 0.1)
    if constraints is None:
        native = meshers.generate(field, bounds=bounds, band=band, **options)
    else:
        if options.pop("order", 1) != 1 or options.pop("threads", 1) != 1:
            raise NotImplementedError(
                "Mapped and clipped TPMS currently require Tet4 and one thread"
            )
        adaptive = options.pop("adapt", False)
        if (adaptive and adaptive_field is None) or options.pop("curved_fit", False):
            raise NotImplementedError(
                "Mapped intersections do not support adaptive or quadratic fitting options"
            )
        options.setdefault("geometry_tolerance", 0.01)
        options.setdefault("max_tetrahedra", 250_000)
        if weld_axis is not None and getattr(meshers, "CLOSED_MAPPED_SEAMS", False):
            options["closed_axes"] = tuple(axis == weld_axis for axis in range(3))
            weld_axis = None
        if adaptive and adaptive_field is not None:
            options.pop("interior_refine_nodes", None)
            options.update(
                adapt=True, inverse_map=inverse_map, physical_bounds=adaptive_field[1]
            )
        native = meshers.generate_intersection(
            constraints,
            bounds=bounds,
            coordinate_map=coordinate_map,
            periodic_transforms=periodic_transforms,
            **options,
        )
        if adaptive and adaptive_field is not None:
            bounds = adaptive_field[1]
    width = native.tetrahedra.shape[1]
    cell_type = {4: pv.CellType.TETRA, 10: pv.CellType.QUADRATIC_TETRA}[width]
    points = native.points.copy()
    lengths = np.asarray(bounds)[1::2] - np.asarray(bounds)[::2]
    shifts = np.diag(lengths)
    if rotation is not None:
        points = rotation.apply(points) + np.asarray(center)
        shifts = rotation.apply(shifts)
    indices = np.arange(len(points))
    faces = native.surface.copy()
    tags = native.boundary_tags.copy()
    if weld_axis is not None:
        pairs = native.periodic_pairs[weld_axis]
        if len(pairs):
            if not np.allclose(
                points[pairs[:, 0]], points[pairs[:, 1]], rtol=0, atol=1e-10
            ):
                raise ValueError("Full-turn seam nodes do not coincide")
            indices[pairs[:, 1]] = pairs[:, 0]
        retained, indices = np.unique(indices, return_inverse=True)
        keep = (tags != 2 * weld_axis + 1) & (tags != 2 * weld_axis + 2)
        faces, tags = faces[keep], tags[keep]
    else:
        retained = np.arange(len(points))
    grid = pv.UnstructuredGrid(
        {cell_type: indices[native.tetrahedra]}, points[retained]
    )

    if native.parameters is not None:
        grid.point_data["meshers_parameters"] = native.parameters[retained].copy()
    if native.constraint_masks is not None:
        grid.point_data["meshers_constraint_masks"] = native.constraint_masks[
            retained
        ].copy()
        grid.field_data["meshers_constraint_names"] = np.asarray(
            native.constraint_names
        )
        edges = indices[native.feature_edges]
        if len(edges):
            grid.field_data["meshers_feature_edges"] = np.unique(
                np.sort(edges, axis=1), axis=0
            )

    if len(faces):
        grid.field_data["meshers_boundary_faces"] = indices[faces]
        grid.field_data["meshers_boundary_tags"] = tags
    for axis, pairs, shift in zip("xyz", native.periodic_pairs, shifts, strict=True):
        if len(pairs) and (weld_axis is None or "xyz".index(axis) != weld_axis):
            grid.field_data[f"periodic_pairs_{axis}"] = indices[pairs]
            if native.periodic_transforms is None:
                transform = np.eye(4)
                transform[:3, 3] = shift
            else:
                placement = np.eye(4)
                if rotation is not None:
                    placement[:3, :3] = rotation.as_matrix()
                    placement[:3, 3] = center
                transform = (
                    placement
                    @ native.periodic_transforms["xyz".index(axis)]
                    @ np.linalg.inv(placement)
                )
            grid.field_data[f"periodic_transform_{axis}"] = transform
            if np.allclose(transform[:3, :3], np.eye(3), atol=1e-12):
                grid.field_data[f"periodic_shift_{axis}"] = transform[:3, 3][
                    None, :
                ].copy()
    return grid
