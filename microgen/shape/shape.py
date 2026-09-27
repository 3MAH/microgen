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
from scipy.spatial.transform import Rotation

from microgen import _meshers

from . import implicit_ops as _ops
from ._types import BoundsType, Field, PeriodType

if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, Vector3DType

class ShellCreationError(Exception):
    """Raised when an OCCT shell cannot be created from a mesh."""


class Shape:
    """Unified shape with optional implicit (F-rep) and CAD representations.

    Every shape has a ``center`` and ``orientation``.  It may also carry an
    implicit scalar field (``_func``) where ``f(x, y, z) < 0`` means *inside*.
    When the implicit field is present, the default :meth:`generate_surface_mesh` and
    :meth:`generate_cad` produce geometry via meshers.  Subclasses
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
        self._mesh_cache = {}

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

    def generate_meshers(
        self,
        bounds=None,
        resolution=50,
        *,
        periodic=(False, False, False),
        **options,
    ):
        """Return a native meshers mesh, including diagnostics and periodic pairs.

        Resolution counts grid points per axis. Fields and bounds use world
        coordinates. Meshers failures propagate without a fallback mesh.
        """
        if self._func is None:
            raise NotImplementedError("No implicit field defined on this shape")
        bounds = self._bounds if bounds is None else bounds
        if bounds is None:
            raise ValueError("Bounds must be provided at construction or meshing")
        # Cache only default backend options. Cancellation and custom evaluators
        # must execute on every call. PyVista conversions copy these arrays.
        key = (tuple(bounds), tuple(np.broadcast_to(resolution, (3,))), tuple(periodic))
        if not options and key in self._mesh_cache:
            return self._mesh_cache[key]
        result = _meshers.generate(
            self._func, bounds, resolution, periodic=periodic, **options
        )
        if not options:
            self._mesh_cache[key] = result
        return result

    def generate_surface_mesh(self, bounds=None, resolution=50, **options):
        """Return the closed meshers solid boundary as PyVista triangles."""
        return _meshers.surface(self.generate_meshers(bounds, resolution, **options))

    def generate_volume_mesh(self, bounds=None, resolution=50, **options):
        """Return a meshers tetrahedral mesh as a PyVista UnstructuredGrid."""
        return _meshers.volume(self.generate_meshers(bounds, resolution, **options))

    def generate_cad(
        self: Shape,
        bounds: BoundsType | None = None,
        resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> CadShape:
        """Generate a CAD shape.

        The default implementation delegates to
        :func:`microgen.cad.shape_to_cad`, which builds a tessellated OCCT BREP
        from the meshers solid boundary. Concrete subclasses
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
