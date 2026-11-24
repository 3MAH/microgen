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
    default_periodic: tuple[bool, bool, bool] = (False, False, False),
    length_scale: float | None = None,
    band=None,
    rotation: Rotation | None = None,
    center=(0.0, 0.0, 0.0),
) -> pv.UnstructuredGrid:
    """Mesh ``field < 0`` (or ``band``) on ``bounds`` with Meshers.

    Tet4/Tet10 output is converted to VTK, placed by ``rotation`` then
    ``center``, and keeps the boundary facets and periodic node pairs in
    ``field_data``.

    :param field: implicit field, negative inside
    :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)`` meshed by Meshers
    :param options: :func:`meshers.generate` options, plus
        ``geometry_tolerance_mode``: ``"relative"`` (default) multiplies
        ``geometry_tolerance`` (Meshers default ``0.01``) by ``length_scale``;
        ``"absolute"`` passes it unchanged, in length units
    :param default_cells: ``cells`` when the options do not set it
    :param default_periodic: ``periodic`` when the options do not set it
    :param length_scale: length of ``"relative"`` tolerances; defaults to
        the smallest extent of ``bounds``
    :param band: ``(lower, upper)``: mesh ``lower < field < upper`` instead
    :param rotation: rotation applied to the points and periodic shifts
    :param center: translation applied after ``rotation``
    :return: tetrahedral mesh
    :raises ImportError: if Meshers is not installed
    :raises ValueError: if ``bounds`` is ``None``, ``options`` holds ``band``
        or ``geometry_tolerance_mode`` is unknown
    """
    if bounds is None:
        raise ValueError("Bounds are required for Meshers volume meshing")
    options = dict(options)
    lengths = np.asarray(bounds, dtype=np.float64)
    lengths = lengths[1::2] - lengths[::2]
    mode = options.pop("geometry_tolerance_mode", "relative")
    if mode == "relative":
        scale = float(lengths.min()) if length_scale is None else length_scale
        options["geometry_tolerance"] = options.get("geometry_tolerance", 0.01) * scale
    elif mode != "absolute":
        raise ValueError(
            f"geometry_tolerance_mode must be 'relative' or 'absolute', not {mode!r}"
        )
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
    options.setdefault("periodic", default_periodic)
    options.setdefault("threads", 1)
    options.setdefault("minimum_quality", 0.1)
    native = meshers.generate(field, bounds=bounds, band=band, **options)
    width = native.tetrahedra.shape[1]
    cell_type = {4: pv.CellType.TETRA, 10: pv.CellType.QUADRATIC_TETRA}[width]
    points, shifts = native.points, np.diag(lengths)
    if rotation is not None:
        points, shifts = rotation.apply(points), rotation.apply(shifts)
    # The sum also detaches the points from the native buffer.
    grid = pv.UnstructuredGrid(
        {cell_type: native.tetrahedra}, points + np.asarray(center)
    )

    if len(native.surface):
        grid.field_data["meshers_boundary_faces"] = native.surface.copy()
        grid.field_data["meshers_boundary_tags"] = native.boundary_tags.copy()
    for axis, pairs, shift in zip("xyz", native.periodic_pairs, shifts, strict=True):
        if len(pairs):
            grid.field_data[f"periodic_pairs_{axis}"] = pairs.copy()
            grid.field_data[f"periodic_shift_{axis}"] = shift[None, :].copy()
    return grid
