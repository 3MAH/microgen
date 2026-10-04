"""Basic Geometry.

====================================================
Basic Geometry (:mod:`microgen.shape.shape`)
====================================================
"""

from __future__ import annotations

import copy as _copy
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.spatial.transform import Rotation

from . import _affine
from . import implicit_ops as _ops
from ._types import BoundsType, Field, PeriodType

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

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
    implicit scalar field (:attr:`field`), in world coordinates, where
    ``f(x, y, z) < 0`` means *inside*.
    When the implicit field is present, the default :meth:`generate_surface_mesh` and
    :meth:`generate_cad` produce geometry via marching cubes.  Subclasses
    (e.g. ``Sphere``, ``Tpms``) override these methods with their own
    implementations.

    Boolean operators (``|``, ``&``, ``-``, ``~``) and smooth boolean
    methods operate on the implicit field and return a new :class:`Shape`.

    :meth:`translate`, :meth:`rotate` and :meth:`scale` follow the PyVista
    convention: ``inplace=False`` (default) returns a new shape,
    ``inplace=True`` mutates ``self``, and both return the shape.  A
    subclass keeps its class and updates its native parameters (``radius``,
    ``dim``, ``cell_size``...) when they can express the transform.  When
    they cannot (e.g. a per-axis scale of a
    :class:`~microgen.shape.sphere.Sphere`), ``inplace=False`` returns a
    generic :class:`Shape` and ``inplace=True`` raises ``ValueError``.

    :param center: center of the shape
    :param orientation: orientation of the shape
    :param field: implicit scalar field ``(x, y, z) -> array``, or ``None``
    :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)`` or ``None``
    :param period: ``(Lx, Ly, Lz)`` if the field is intrinsically periodic
        (``field(p + L) == field(p)`` along each axis), or ``None``.
        Set by ``Tpms`` and ``Spinodoid`` from ``cell_size * repeat_cell``.
    """

    def __init__(
        self: Shape,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType | Rotation = (0, 0, 0),
        field: Field | None = None,
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
        self._field = field
        self._bounds = bounds
        self._period: PeriodType | None = period
        # Cache of sampled structured grids keyed on (bounds, resolution).
        # Shared between generate_surface_mesh and generate_volume_mesh so
        # users calling both on the same instance only pay one N^3 field
        # evaluation. Cleared by every in-place transform.
        self._grid_cache: dict[tuple[BoundsType, int], pv.StructuredGrid] = {}

    # ------------------------------------------------------------------
    # Public read-only accessors
    # ------------------------------------------------------------------

    @property
    def center(self: Shape) -> Vector3DType:
        """Geometric center (read-only; changed by the transforms).

        Subclasses with a native renderer (``Sphere``, ``Tpms``, …) read
        this in their ``generate_*`` overrides and bake it into the field.
        Use :meth:`translate`, :meth:`rotate` or :meth:`scale` to move a
        shape: they keep the field, the native parameters and the
        renderers coherent.
        """
        return self._center

    @property
    def orientation(self: Shape) -> Rotation:
        """Rotation applied by subclasses' renderers (read-only; see :meth:`rotate`)."""
        return self._orientation

    @property
    def field(self: Shape) -> Field | None:
        """The implicit scalar field, or ``None``."""
        return self._field

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

    def require_field(self: Shape) -> Field:
        """Return the implicit field or raise ``ValueError`` if it is not set."""
        if self._field is None:
            err_msg = "No implicit scalar field defined on this shape"
            raise ValueError(err_msg)
        return self._field

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

        Coordinates are world coordinates: subclasses bake ``center`` and
        ``orientation`` into the field, and :meth:`translate`,
        :meth:`rotate` and :meth:`scale` compose with it.

        :param x: x coordinates
        :param y: y coordinates
        :param z: z coordinates
        :return: scalar field values (negative = inside)
        """
        return self.require_field()(x, y, z)

    # ------------------------------------------------------------------
    # Mesh generation (defaults use the implicit field)
    # ------------------------------------------------------------------

    def _sample_implicit_grid(
        self: Shape,
        bounds: BoundsType | None,
        resolution: int,
        caller: str,
    ) -> pv.StructuredGrid:
        """Build a structured grid over ``bounds`` and sample the field onto it.

        Shared by :meth:`generate_surface_mesh` and :meth:`generate_volume_mesh`,
        with a per-instance ``(bounds, resolution)`` cache so consecutive calls
        on the same shape only pay one N^3 field evaluation. Raises
        ``NotImplementedError`` (with a caller-specific message) when the field
        is unset, and ``ValueError`` when bounds can't be resolved.
        """
        if self._field is None:
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
        with non-zero ``center`` / ``orientation`` bake those into the field
        during construction).

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
    # Copy and transforms — PyVista convention: imperative verb +
    # ``inplace=False``.
    #
    # Subclasses with native parameters (``center``, ``orientation``,
    # ``radius``, ``dim``, ``cell_size``…) are the source of truth: a
    # transform they can express updates those parameters and rebuilds the
    # field, so field, CAD and meshes stay coherent.  A transform they
    # cannot express returns a generic :class:`Shape` (``inplace=False``)
    # or raises (``inplace=True``, the class cannot change in place).
    # ------------------------------------------------------------------

    #: ``True`` on subclasses whose field is rebuilt from native parameters
    #: by :meth:`_rebuild_field`.
    _native_params: bool = False

    def copy(self: Shape) -> Shape:
        """Return a shallow copy that keeps the subclass and its parameters.

        The field callable and immutable parameters are shared; the sampled
        grid cache and any mutable state listed by the subclass are not, so
        an ``inplace=True`` transform on the copy never touches ``self``.
        """
        new = _copy.copy(self)
        new._grid_cache = {}  # noqa: SLF001
        new._post_copy()  # noqa: SLF001
        return new

    def _post_copy(self: Shape) -> None:
        """Detach mutable state a shallow copy would share (subclass hook)."""

    def _rebuild_field(self: Shape) -> None:
        """Rebuild ``_field`` / ``_bounds`` / ``_period`` from native parameters."""
        self._setup_frep_field()  # type: ignore[attr-defined]

    def _scale_params(self: Shape, factors: npt.NDArray[np.float64]) -> bool:
        """Rescale native size parameters by world ``factors`` (subclass hook).

        Return ``False``, leaving ``self`` untouched, when the scale cannot be
        expressed in the native parameters.
        """
        return False

    def _local_scale_factors(
        self: Shape,
        factors: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64] | None:
        """World scale factors seen along the local axes, or ``None``.

        Defined when the scale is uniform or when ``orientation`` is a
        signed permutation (the local axes lie on world axes).
        """
        if np.all(factors == factors[0]):
            return factors.copy()
        perm = _affine.signed_permutation(self._orientation.as_matrix())
        if perm is None:
            return None
        return np.abs(perm).T @ factors

    def _posed(
        self: Shape,
        field: Field,
        bounds: BoundsType | None,
    ) -> tuple[Field, BoundsType | None]:
        """Map a local-frame field and bounds to world through center and orientation."""
        center = np.asarray(self._center, dtype=np.float64)
        matrix = self._orientation.as_matrix()
        if not center.any() and np.array_equal(matrix, np.eye(3)):
            return field, bounds
        zero = np.zeros(3)
        return (
            _affine.transform_field(field, matrix, zero, center),
            None
            if bounds is None
            else _affine.transform_bounds(bounds, matrix, zero, center),
        )

    def _set_local_frep(
        self: Shape,
        field: Field,
        bounds: BoundsType | None,
        period: PeriodType | None,
    ) -> None:
        """Store a local-frame field and pose it in world (``Tpms``, ``Spinodoid``).

        Classes whose renderers work in a local frame (a cached grid rotated
        and translated at the end) keep the local field for those renderers
        and expose the world field, so ``field`` agrees with the meshes and
        the CAD for any ``center`` and ``orientation``.
        """
        self._local_field = field
        self._local_bounds = bounds
        self._local_period = period
        self._pose_local_frep()

    def _pose_local_frep(self: Shape) -> None:
        """Set ``_field`` / ``_bounds`` / ``_period`` from the local-frame ones."""
        self._field, self._bounds = self._posed(self._local_field, self._local_bounds)
        self._period = _affine.rotate_period(self._local_period, self._orientation)

    def translate(
        self: Shape,
        offset: Sequence[float],
        *,
        inplace: bool = False,
    ) -> Shape:
        """Translate the shape by *offset* (PyVista convention).

        ``evaluate(p) == old.evaluate(p - offset)``; ``center`` and
        ``bounds`` shift by *offset*, ``orientation`` and ``period`` are
        kept.

        :param offset: ``(dx, dy, dz)`` world-space shift
        :param inplace: mutate ``self`` when ``True``; otherwise (default)
            return a new shape of the same class
        :return: the translated shape
        """
        shift = np.asarray(offset, dtype=np.float64).reshape(3)
        zero = np.zeros(3)

        def _update(target: Shape) -> bool:
            target._center = _affine.transform_point(  # noqa: SLF001
                target._center,  # noqa: SLF001
                np.eye(3),
                zero,
                shift,
            )
            return True

        return self._transform(
            np.eye(3),
            zero,
            shift,
            rotation=None,
            period=self._period,
            update_params=_update,
            inplace=inplace,
        )

    def rotate(
        self: Shape,
        rotation: Rotation | npt.ArrayLike,
        point: Sequence[float] | None = None,
        *,
        inplace: bool = False,
    ) -> Shape:
        """Rotate the shape about *point* (PyVista convention).

        ``center`` rotates about *point*, ``orientation`` composes on the
        left with *rotation* and ``bounds`` becomes the AABB of the rotated
        AABB.  A rotation within ``1e-12`` of a signed permutation (quarter
        turns about the axes) is snapped to it and permutes ``period``;
        any other rotation sets ``period`` to ``None``.

        :param rotation: a SciPy ``Rotation`` or a proper 3x3 rotation matrix
        :param point: pivot; defaults to the world origin
        :param inplace: mutate ``self`` when ``True``; otherwise (default)
            return a new shape of the same class
        :return: the rotated shape
        """
        rot = _affine.as_rotation(rotation)
        matrix = rot.as_matrix()
        pivot = _affine.as_point(point)
        zero = np.zeros(3)

        def _update(target: Shape) -> bool:
            target._center = _affine.transform_point(  # noqa: SLF001
                target._center,  # noqa: SLF001
                matrix,
                pivot,
                zero,
            )
            target._orientation = rot * target._orientation  # noqa: SLF001
            return True

        return self._transform(
            matrix,
            pivot,
            zero,
            rotation=rot,
            period=_affine.rotate_period(self._period, rot),
            update_params=_update,
            inplace=inplace,
        )

    def scale(
        self: Shape,
        factor: float | Sequence[float],
        point: Sequence[float] | None = None,
        *,
        inplace: bool = False,
    ) -> Shape:
        """Scale the shape about *point* (PyVista convention).

        A uniform factor keeps a signed-distance field exact
        (``s * f(x / s)``).  A per-axis factor multiplies the composed field
        by ``min(s)``: the solid is exact but the field is only a lower
        bound of the distance.  ``period`` scales with the factors.

        :param factor: uniform factor or ``(sx, sy, sz)``, strictly positive
        :param point: pivot; defaults to the world origin
        :param inplace: mutate ``self`` when ``True``; otherwise (default)
            return a new shape (of the same class when its native parameters
            can express the scale, a generic :class:`Shape` otherwise)
        :return: the scaled shape
        :raises ValueError: if a factor is not strictly positive, or if
            ``inplace=True`` and the native parameters cannot express the
            scale
        """
        factors = _affine.scale_factors(factor)
        matrix = np.diag(factors)
        pivot = _affine.as_point(point)
        zero = np.zeros(3)

        def _update(target: Shape) -> bool:
            if not target._scale_params(factors):  # noqa: SLF001
                return False
            target._center = _affine.transform_point(  # noqa: SLF001
                target._center,  # noqa: SLF001
                matrix,
                pivot,
                zero,
            )
            return True

        return self._transform(
            matrix,
            pivot,
            zero,
            rotation=None,
            period=_affine.scale_period(self._period, factors),
            update_params=_update,
            inplace=inplace,
            multiplier=float(factors.min()),
        )

    def _transform(  # noqa: PLR0913
        self: Shape,
        matrix: npt.NDArray[np.float64],
        point: npt.NDArray[np.float64],
        shift: npt.NDArray[np.float64],
        *,
        rotation: Rotation | None,
        period: PeriodType | None,
        update_params: Callable[[Shape], bool],
        inplace: bool,
        multiplier: float = 1.0,
    ) -> Shape:
        """Apply ``x -> matrix (x - point) + point + shift`` (shared worker)."""
        field = self.require_field()
        if self._native_params:
            target = self if inplace else self.copy()
            if update_params(target):
                target._rebuild_field()  # noqa: SLF001
                target._grid_cache.clear()  # noqa: SLF001
                return target
            if inplace:
                err_msg = (
                    f"{type(self).__name__} cannot express this transform in "
                    "its native parameters; use inplace=False to get a "
                    "generic Shape"
                )
                raise ValueError(err_msg)

        new_field = _affine.transform_field(field, matrix, point, shift, multiplier)
        new_bounds = (
            None
            if self._bounds is None
            else _affine.transform_bounds(self._bounds, matrix, point, shift)
        )
        new_center = _affine.transform_point(self._center, matrix, point, shift)
        new_orientation = (
            self._orientation if rotation is None else rotation * self._orientation
        )
        if self._native_params:
            return Shape(
                center=new_center,
                orientation=new_orientation,
                field=new_field,
                bounds=new_bounds,
                period=period,
            )
        target = self if inplace else self.copy()
        target._field = new_field  # noqa: SLF001
        target._bounds = new_bounds  # noqa: SLF001
        target._center = new_center  # noqa: SLF001
        target._orientation = new_orientation  # noqa: SLF001
        target._period = period  # noqa: SLF001
        target._grid_cache.clear()  # noqa: SLF001
        return target
