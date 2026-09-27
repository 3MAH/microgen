"""Basic Geometry.

====================================================
Basic Geometry (:mod:`microgen.shape.shape`)
====================================================
"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.spatial.transform import Rotation

from . import implicit_ops as _ops
from ._types import BoundsType, Field, PeriodType

if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, Vector3DType

# Scalar-field name on the StructuredGrid used by the default mesh generators.
_IMPLICIT_SCALAR = "implicit"


def _pad_grid_with_outside_halo(
    grid: pv.StructuredGrid,
    scalar: str,
) -> pv.StructuredGrid:
    """Wrap ``grid`` in a single-cell halo with a large positive scalar value.

    Used by :meth:`Shape.generate_surface_mesh` so that ``contour(0.0)`` closes
    the iso-surface at the original bbox: marching cubes finds zeros between
    interior negative nodes and the halo's "outside" nodes, producing cap
    triangles AT the original bounds rather than leaving them open.
    """
    nx, ny, nz = grid.dimensions
    pts = np.asarray(grid.points).reshape((nx, ny, nz, 3), order="F")
    dx = float(pts[1, 0, 0, 0] - pts[0, 0, 0, 0])
    dy = float(pts[0, 1, 0, 1] - pts[0, 0, 0, 1])
    dz = float(pts[0, 0, 1, 2] - pts[0, 0, 0, 2])
    xi = np.concatenate(
        [[pts[0, 0, 0, 0] - dx], pts[:, 0, 0, 0], [pts[-1, 0, 0, 0] + dx]]
    )
    yi = np.concatenate(
        [[pts[0, 0, 0, 1] - dy], pts[0, :, 0, 1], [pts[0, -1, 0, 1] + dy]]
    )
    zi = np.concatenate(
        [[pts[0, 0, 0, 2] - dz], pts[0, 0, :, 2], [pts[0, 0, -1, 2] + dz]]
    )
    x, y, z = np.meshgrid(xi, yi, zi, indexing="ij")

    field = np.asarray(grid[scalar]).reshape((nx, ny, nz), order="F")
    pad_val = max(float(np.nanmax(field)) + 1.0, 1e6)
    padded = np.full((nx + 2, ny + 2, nz + 2), pad_val, dtype=field.dtype)
    padded[1:-1, 1:-1, 1:-1] = field

    out = pv.StructuredGrid(x, y, z)
    out[scalar] = padded.ravel(order="F")
    return out


class ShellCreationError(Exception):
    """Raised when an OCCT shell cannot be created from a mesh."""


class Shape:
    """Unified shape with optional implicit (F-rep) and CAD representations.

    Every shape has a ``center`` and ``orientation``.  It may also carry an
    implicit scalar field (``_func``) where ``f(x, y, z) < 0`` means *inside*.
    When the implicit field is present, the default :meth:`generate_surface_mesh` and
    :meth:`generate_cad` produce geometry via marching cubes.  Subclasses
    (e.g. ``Sphere``, ``Tpms``) override these methods with their own
    implementations.

    Boolean operators (``|``, ``&``, ``-``, ``~``) and smooth boolean
    methods operate on the implicit field and return a new :class:`Shape`.

    :param center: center of the shape
    :param orientation: orientation of the shape
    :param func: implicit scalar field ``(x, y, z) -> array``, or ``None``
    :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)`` or ``None``
    :param period: ``(Lx, Ly, Lz)`` if the field is intrinsically periodic
        (``func(p + L) == func(p)`` along each axis), or ``None``.
        Set by ``Tpms`` and ``Spinodoid`` from ``cell_size * repeat_cell``.
    """

    def __init__(
        self: Shape,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType | Rotation = (0, 0, 0),
        func: Field | None = None,
        bounds: BoundsType | None = None,
        period: PeriodType | None = None,
    ) -> None:
        """Initialize the shape."""
        self._center = center
        self._orientation = (
            orientation
            if isinstance(orientation, Rotation)
            else Rotation.from_euler("ZXZ", orientation, degrees=True)
        )
        self._func = func
        self._bounds = bounds
        self._period: PeriodType | None = period
        # Cache of sampled structured grids keyed on (bounds, resolution).
        # Shared between generate_surface_mesh and generate_volume_mesh so
        # users calling both on the same instance only pay one N^3 field
        # evaluation. Cleared lazily — Shape is immutable post-construction
        # (center/orientation are read-only properties, _func is fixed).
        self._grid_cache: dict[tuple[BoundsType, int], pv.StructuredGrid] = {}

    # ------------------------------------------------------------------
    # Public read-only accessors
    # ------------------------------------------------------------------

    @property
    def center(self: Shape) -> Vector3DType:
        """Geometric center (set at construction, immutable).

        Subclasses with a native renderer (``Sphere``, ``Tpms``, …) read
        this in their ``generate_*`` overrides. For an implicit-only
        :class:`Shape` (built via :func:`microgen.shape.implicit_ops.from_field`
        or by a boolean composition), use :meth:`translate` to shift the
        field — assigning to ``.center`` after construction has no effect
        on a bare ``Shape``, so the attribute is read-only at the base.
        """
        return self._center

    @property
    def orientation(self: Shape) -> Rotation:
        """Rotation applied by subclasses' renderers (set at construction, immutable).

        Implicit-only shapes should compose with :meth:`rotate` instead
        of mutating this attribute.
        """
        return self._orientation

    @property
    def func(self: Shape) -> Field | None:
        """The implicit scalar field, or ``None``."""
        return self._func

    @property
    def bounds(self: Shape) -> BoundsType | None:
        """The bounding box ``(xmin, xmax, ymin, ymax, zmin, zmax)``, or ``None``."""
        return self._bounds

    @property
    def period(self: Shape) -> PeriodType | None:
        """The intrinsic period ``(Lx, Ly, Lz)`` if the field is periodic, else ``None``.

        When non-``None``, ``self.evaluate(x + Lx, y, z) == self.evaluate(x, y, z)``
        (and analogously for y, z) — i.e. periodicity is a data-structure
        invariant of the field, not a runtime flag.  ``Tpms`` and
        ``Spinodoid`` set this from ``cell_size * repeat_cell``.
        """
        return self._period

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
        """Evaluate the implicit scalar field at the given coordinates.

        Coordinates are in the **field's local frame** — ``center`` and
        ``orientation`` are NOT applied here (they only affect mesh output
        in :meth:`generate_surface_mesh`).  Use :meth:`translate` / :meth:`rotate`
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

    def _sample_implicit_grid(
        self: Shape,
        bounds: BoundsType | None,
        resolution: int,
        caller: str,
    ) -> pv.StructuredGrid:
        """Build a structured grid over ``bounds`` and sample ``_func`` onto it.

        Shared by :meth:`generate_surface_mesh` and :meth:`generate_volume_mesh`,
        with a per-instance ``(bounds, resolution)`` cache so consecutive calls
        on the same shape only pay one N^3 field evaluation. Raises
        ``NotImplementedError`` (with a caller-specific message) when ``_func``
        is unset, and ``ValueError`` when bounds can't be resolved.
        """
        if self._func is None:
            err_msg = f"No implicit field defined — subclasses must override {caller}()"
            raise NotImplementedError(err_msg)

        bounds = bounds or self._bounds
        if bounds is None:
            err_msg = f"Bounds must be provided either at construction or in {caller}()"
            raise ValueError(err_msg)

        cache_key = (tuple(bounds), resolution)
        cached = self._grid_cache.get(cache_key)
        if cached is not None:
            return cached

        xmin, xmax, ymin, ymax, zmin, zmax = bounds
        xi = np.linspace(xmin, xmax, resolution)
        yi = np.linspace(ymin, ymax, resolution)
        zi = np.linspace(zmin, zmax, resolution)
        x, y, z = np.meshgrid(xi, yi, zi, indexing="ij")

        grid = pv.StructuredGrid(x, y, z)
        grid[_IMPLICIT_SCALAR] = self.evaluate(
            x.ravel(order="F"),
            y.ravel(order="F"),
            z.ravel(order="F"),
        )
        self._grid_cache[cache_key] = grid
        return grid

    def generate_surface_mesh(
        self: Shape,
        bounds: BoundsType | None = None,
        resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate a surface VTK mesh of the shape.

        The default implementation runs marching cubes (``f < 0``) on the
        cached implicit grid wrapped in a single-cell halo of "outside"
        values, so the iso-surface naturally closes at the bbox: where the
        volume reaches the bounds, cap faces are produced AT the bbox.
        Subclasses with a native renderer (``Sphere``, ``Box``, ``Tpms``, …)
        override this.

        The implicit field is expected to be in world coordinates (subclasses
        with non-zero ``center`` / ``orientation`` should bake those into
        ``_func`` during construction).

        The sampled structured grid is cached per ``(bounds, resolution)``
        on the instance, shared with :meth:`generate_volume_mesh`. The cache
        is unbounded — calling this method with many distinct ``resolution``
        values on the same instance retains every sampled grid until the
        instance is GC'd.

        :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)``
        :param resolution: number of grid points per axis
        :return: triangulated surface mesh
        """
        grid = self._sample_implicit_grid(bounds, resolution, "generate_surface_mesh")
        padded = _pad_grid_with_outside_halo(grid, _IMPLICIT_SCALAR)
        polydata = padded.contour(isosurfaces=[0.0], scalars=_IMPLICIT_SCALAR)
        if polydata.n_cells == 0:
            return pv.PolyData()
        return polydata.clean().triangulate()

    def generate_volume_mesh(
        self: Shape,
        bounds: BoundsType | None = None,
        resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.UnstructuredGrid:
        """Generate a volumetric VTK mesh of the shape's interior.

        Default implementation samples the implicit field on a structured
        grid over ``bounds`` and keeps cells where ``f < 0``. Subclasses
        with a native volumetric representation (``Tpms``, ``Spinodoid``)
        override this with their cached-grid path.

        The implicit field is expected to be in world coordinates (same
        contract as :meth:`generate_surface_mesh`).

        :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)``
        :param resolution: number of grid points per axis
        :return: clipped ``pv.UnstructuredGrid`` covering the shape's interior
        """
        grid = self._sample_implicit_grid(bounds, resolution, "generate_volume_mesh")
        return grid.clip_scalar(scalars=_IMPLICIT_SCALAR, value=0.0, invert=True)

    def generate_cad(
        self: Shape,
        bounds: BoundsType | None = None,
        resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> CadShape:
        """Generate a CAD shape.

        The default implementation delegates to
        :func:`microgen.cad.shape_to_cad`, which builds a tessellated OCCT BREP
        from the implicit-field marching-cubes mesh.  Concrete subclasses
        with native primitive paths (``Box``, ``Sphere``, ``Cylinder``,
        ``Capsule``, ``Ellipsoid``, ``Tpms``, ``Spinodoid`` …) override
        this method with native OCCT construction.

        Requires the optional ``[cad]`` install extra (``cadquery-ocp-novtk``).

        :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)``
        :param resolution: number of grid points per axis
        :return: :class:`microgen.cad.CadShape` wrapping an OCCT ``TopoDS_Shell``
        """
        from microgen.cad import shape_to_cad  # noqa: PLC0415

        return shape_to_cad(self, bounds=bounds, resolution=resolution)

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
        """Return a new shape translated by *offset*.

        The returned :class:`Shape` has its ``center`` shifted by *offset*
        and its ``bounds`` updated; ``orientation`` is preserved. The
        implicit field is composed so ``evaluate(p) == old.evaluate(p - offset)``.
        """
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
        cx, cy, cz = self._center
        return Shape(
            func=lambda x, y, z, _f=f, _dx=dx, _dy=dy, _dz=dz: _f(
                x - _dx,
                y - _dy,
                z - _dz,
            ),
            bounds=new_bounds,
            center=(cx + dx, cy + dy, cz + dz),
            orientation=self._orientation,
        )

    def rotate(
        self: Shape,
        angles: tuple[float, float, float],
        convention: str = "ZXZ",
    ) -> Shape:
        """Return a new shape rotated by Euler *angles* (degrees).

        The rotation is applied **about the world origin**. The returned
        shape's ``center`` is the rotated original center, ``orientation``
        composes left with the rotation, and ``bounds`` is the AABB of
        the rotated original AABB.
        """
        f = self.require_func()
        rot = Rotation.from_euler(convention, angles, degrees=True)
        rot_matrix = rot.as_matrix()
        inv_matrix = rot.inv().as_matrix()
        new_bounds = None
        if self._bounds is not None:
            b = self._bounds
            corners = np.array(
                list(itertools.product(b[0:2], b[2:4], b[4:6])),
            )
            rotated = (rot_matrix @ corners.T).T
            new_bounds = (
                float(rotated[:, 0].min()),
                float(rotated[:, 0].max()),
                float(rotated[:, 1].min()),
                float(rotated[:, 1].max()),
                float(rotated[:, 2].min()),
                float(rotated[:, 2].max()),
            )
        rotated_center = rot_matrix @ np.asarray(self._center, dtype=np.float64)
        return Shape(
            func=lambda x, y, z, _f=f, _m=inv_matrix: _f(
                *(_m @ np.array([x, y, z])),
            ),
            bounds=new_bounds,
            center=tuple(rotated_center.tolist()),
            orientation=rot * self._orientation,
        )

    def scale(self: Shape, factor: float) -> Shape:
        """Return a new shape uniformly scaled by *factor* about the world origin.

        ``center`` is scaled by the same factor; ``orientation`` is
        preserved; ``bounds`` is scaled (with axis-pair swap for
        negative factors).
        """
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
                new_bounds = (
                    new_bounds[1],
                    new_bounds[0],
                    new_bounds[3],
                    new_bounds[2],
                    new_bounds[5],
                    new_bounds[4],
                )
        cx, cy, cz = self._center
        return Shape(
            func=lambda x, y, z, _f=f, _s=factor: _f(x / _s, y / _s, z / _s) * _s,
            bounds=new_bounds,
            center=(cx * factor, cy * factor, cz * factor),
            orientation=self._orientation,
        )
