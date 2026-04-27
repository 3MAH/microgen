"""
Basic Geometry.

====================================================
Basic Geometry (:mod:`microgen.shape.shape`)
====================================================
"""

from __future__ import annotations

import itertools
from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.spatial.transform import Rotation

from microgen.operations import rotate as rotate_mesh

from . import implicit_ops as _ops

if TYPE_CHECKING:
    import cadquery as cq

    from microgen.shape import KwargsGenerateType, Vector3DType

Field = Callable[
    [npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]],
    npt.NDArray[np.float64],
]

BoundsType = tuple[float, float, float, float, float, float]


class ShellCreationError(Exception):
    """Raised when a CadQuery shell cannot be created from a mesh."""


class Shape:
    """
    Unified shape with optional implicit (F-rep) and CAD representations.

    Every shape has a ``center`` and ``orientation``.  It may also carry an
    implicit scalar field (``_func``) where ``f(x, y, z) < 0`` means *inside*.
    When the implicit field is present, the default :meth:`generate_vtk` and
    :meth:`generate` produce geometry via marching cubes.  Subclasses
    (e.g. ``Sphere``, ``Tpms``) override these methods with their own
    implementations.

    Boolean operators (``|``, ``&``, ``-``, ``~``) and smooth boolean
    methods operate on the implicit field and return a new :class:`Shape`.

    :param center: center of the shape
    :param orientation: orientation of the shape
    :param func: implicit scalar field ``(x, y, z) -> array``, or ``None``
    :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)`` or ``None``
    """

    def __init__(
        self: Shape,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType | Rotation = (0, 0, 0),
        func: Field | None = None,
        bounds: BoundsType | None = None,
    ) -> None:
        """Initialize the shape."""
        self.center = center
        self.orientation = (
            orientation
            if isinstance(orientation, Rotation)
            else Rotation.from_euler("ZXZ", orientation, degrees=True)
        )
        self._func = func
        self._bounds = bounds

    # ------------------------------------------------------------------
    # Public read-only accessors for implicit field
    # ------------------------------------------------------------------

    @property
    def func(self: Shape) -> Field | None:
        """The implicit scalar field, or ``None``."""
        return self._func

    @property
    def bounds(self: Shape) -> BoundsType | None:
        """The bounding box ``(xmin, xmax, ymin, ymax, zmin, zmax)``, or ``None``."""
        return self._bounds

    def require_func(self: Shape) -> Field:
        """Return ``_func`` or raise if not set."""
        if self._func is None:
            err_msg = "No implicit scalar field defined on this shape"
            raise ValueError(err_msg)
        return self._func

    # ------------------------------------------------------------------
    # Implicit field evaluation
    # ------------------------------------------------------------------

    def evaluate(
        self: Shape,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """
        Evaluate the implicit scalar field at the given coordinates.

        Coordinates are in the **field's local frame** — ``center`` and
        ``orientation`` are NOT applied here (they only affect mesh output
        in :meth:`generate_vtk`).  Use :meth:`translate` / :meth:`rotate`
        to bake transforms into the field itself.

        :param x: x coordinates
        :param y: y coordinates
        :param z: z coordinates
        :return: scalar field values (negative = inside)
        """
        return self.require_func()(x, y, z)

    # ------------------------------------------------------------------
    # Mesh generation (defaults use the implicit field)
    # ------------------------------------------------------------------

    def generate_vtk(
        self: Shape,
        bounds: BoundsType | None = None,
        resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """
        Generate a VTK mesh of the shape.

        The default implementation meshes the implicit field via marching cubes
        (``f < 0`` convention).  Subclasses override this with their own
        geometry generation.

        :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)``
        :param resolution: number of grid points per axis
        :return: triangulated surface mesh
        """
        if self._func is None:
            err_msg = (
                "No implicit field defined — subclasses must override generate_vtk()"
            )
            raise NotImplementedError(err_msg)

        bounds = bounds or self._bounds
        if bounds is None:
            err_msg = (
                "Bounds must be provided either at construction or in generate_vtk()"
            )
            raise ValueError(err_msg)

        xmin, xmax, ymin, ymax, zmin, zmax = bounds
        xi = np.linspace(xmin, xmax, resolution)
        yi = np.linspace(ymin, ymax, resolution)
        zi = np.linspace(zmin, zmax, resolution)
        x, y, z = np.meshgrid(xi, yi, zi, indexing="ij")

        grid = pv.StructuredGrid(x, y, z)
        field = self.evaluate(
            x.ravel(order="F"),
            y.ravel(order="F"),
            z.ravel(order="F"),
        )
        grid["implicit"] = field

        polydata = grid.contour(isosurfaces=[0.0], scalars="implicit")
        if polydata.n_cells == 0:
            return pv.PolyData()

        polydata = polydata.clean().triangulate()

        polydata = rotate_mesh(polydata, center=(0, 0, 0), rotation=self.orientation)
        return polydata.translate(xyz=self.center)

    def generate(
        self: Shape,
        bounds: BoundsType | None = None,
        resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> cq.Shape:
        """
        Generate a CAD shape.

        The default implementation builds a CadQuery shape from the
        implicit-field VTK mesh.  Subclasses override this with native
        CAD construction.

        :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)``
        :param resolution: number of grid points per axis
        :return: CadQuery Shape
        """
        if self._func is None:
            err_msg = "No implicit field defined — subclasses must override generate()"
            raise NotImplementedError(err_msg)

        import cadquery as cq

        mesh = self.generate_vtk(bounds=bounds, resolution=resolution)
        if mesh.n_cells == 0:
            err_msg = "Generated mesh is empty — check bounds and field function"
            raise ValueError(err_msg)

        if not mesh.is_all_triangles:
            mesh.triangulate(inplace=True)
        triangles = mesh.faces.reshape(-1, 4)[:, 1:]
        triangles = np.c_[triangles, triangles[:, 0]]

        faces = []
        for tri in triangles:
            lines = [
                cq.Edge.makeLine(
                    cq.Vector(*mesh.points[start]),
                    cq.Vector(*mesh.points[end]),
                )
                for start, end in itertools.pairwise(tri)
            ]
            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))

        try:
            shell = cq.Shell.makeShell(faces)
        except ValueError as err:
            err_msg = (
                "Failed to create the shell, "
                "try to increase the resolution or adjust bounds."
            )
            raise ShellCreationError(err_msg) from err

        return cq.Shape(shell.wrapped)

    def generateVtk(self: Shape, **kwargs: KwargsGenerateType) -> pv.PolyData:  # noqa: N802
        """Deprecated. Use :meth:`generate_vtk` instead."""
        return self.generate_vtk(**kwargs)

    # ------------------------------------------------------------------
    # Boolean operators (on implicit field)
    # ------------------------------------------------------------------

    def __or__(self: Shape, other: Shape) -> Shape:
        """Union (``a | b``): inside where either field is negative."""
        return _ops.union(self, other)

    def __and__(self: Shape, other: Shape) -> Shape:
        """Intersection (``a & b``): inside where both fields are negative."""
        return _ops.intersection(self, other)

    def __sub__(self: Shape, other: Shape) -> Shape:
        """Difference (``a - b``): inside *a* but not *b*."""
        return _ops.difference(self, other)

    def __invert__(self: Shape) -> Shape:
        """Complement (``~a``): negate the field."""
        return _ops.complement(self)

    # ------------------------------------------------------------------
    # Smooth booleans
    # ------------------------------------------------------------------

    def smooth_union(self: Shape, other: Shape, k: float) -> Shape:
        """Smooth union with blending radius *k*."""
        return _ops.smooth_union(self, other, k)

    def smooth_intersection(self: Shape, other: Shape, k: float) -> Shape:
        """Smooth intersection with blending radius *k*."""
        return _ops.smooth_intersection(self, other, k)

    def smooth_difference(self: Shape, other: Shape, k: float) -> Shape:
        """Smooth difference with blending radius *k*."""
        return _ops.smooth_difference(self, other, k)

    # ------------------------------------------------------------------
    # Implicit field transforms
    # ------------------------------------------------------------------

    def translate(self: Shape, offset: tuple[float, float, float]) -> Shape:
        """Return a new shape translated by *offset* (implicit field)."""
        f = self.require_func()
        dx, dy, dz = offset
        new_bounds = None
        if self._bounds is not None:
            b = self._bounds
            new_bounds = (
                b[0] + dx,
                b[1] + dx,
                b[2] + dy,
                b[3] + dy,
                b[4] + dz,
                b[5] + dz,
            )
        return Shape(
            func=lambda x, y, z, _f=f, _dx=dx, _dy=dy, _dz=dz: _f(
                x - _dx,
                y - _dy,
                z - _dz,
            ),
            bounds=new_bounds,
        )

    def rotate(
        self: Shape,
        angles: tuple[float, float, float],
        convention: str = "ZXZ",
    ) -> Shape:
        """Return a new shape rotated by Euler *angles* (degrees, implicit field)."""
        f = self.require_func()
        rot = Rotation.from_euler(convention, angles, degrees=True)
        inv_matrix = rot.inv().as_matrix()
        # Recompute AABB by rotating the 8 corners of the original box
        new_bounds = None
        if self._bounds is not None:
            b = self._bounds
            corners = np.array(
                list(itertools.product(b[0:2], b[2:4], b[4:6])),
            )
            rotated = (rot.as_matrix() @ corners.T).T
            new_bounds = (
                float(rotated[:, 0].min()),
                float(rotated[:, 0].max()),
                float(rotated[:, 1].min()),
                float(rotated[:, 1].max()),
                float(rotated[:, 2].min()),
                float(rotated[:, 2].max()),
            )
        return Shape(
            func=lambda x, y, z, _f=f, _m=inv_matrix: _f(
                *(_m @ np.array([x, y, z])),
            ),
            bounds=new_bounds,
        )

    def scale(self: Shape, factor: float) -> Shape:
        """Return a new shape uniformly scaled by *factor* (implicit field)."""
        f = self.require_func()
        new_bounds = None
        if self._bounds is not None:
            b = self._bounds
            new_bounds = (
                b[0] * factor,
                b[1] * factor,
                b[2] * factor,
                b[3] * factor,
                b[4] * factor,
                b[5] * factor,
            )
            if factor < 0:
                # Negative factor inverts min/max — swap each axis pair
                new_bounds = (
                    new_bounds[1],
                    new_bounds[0],
                    new_bounds[3],
                    new_bounds[2],
                    new_bounds[5],
                    new_bounds[4],
                )
        return Shape(
            func=lambda x, y, z, _f=f, _s=factor: _f(x / _s, y / _s, z / _s) * _s,
            bounds=new_bounds,
        )
