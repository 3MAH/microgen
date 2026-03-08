"""F-rep Implicit Modeling.

====================================================
Implicit Shape (:mod:`microgen.shape.implicit_shape`)
====================================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.spatial.transform import Rotation

from microgen.operations import rotate as rotate_mesh

from .shape import Shape

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, Vector3DType

Field = Callable[
    [npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]],
    npt.NDArray[np.float64],
]

BoundsType = tuple[float, float, float, float, float, float]


def _smooth_min(
    a: npt.NDArray[np.float64],
    b: npt.NDArray[np.float64],
    k: float,
) -> npt.NDArray[np.float64]:
    """Smooth minimum (Inigo Quilez cubic polynomial)."""
    if k <= 0:
        return np.minimum(a, b)
    h = np.maximum(k - np.abs(a - b), 0.0) / k
    return np.minimum(a, b) - h * h * h * k / 6.0


def _smooth_max(
    a: npt.NDArray[np.float64],
    b: npt.NDArray[np.float64],
    k: float,
) -> npt.NDArray[np.float64]:
    """Smooth maximum."""
    return -_smooth_min(-a, -b, k)


def _merge_bounds(
    a: BoundsType | None,
    b: BoundsType | None,
    mode: str = "union",
) -> BoundsType | None:
    """Merge two bounding boxes."""
    if a is None and b is None:
        return None
    if a is None:
        return b
    if b is None:
        return a
    if mode == "union":
        return (
            min(a[0], b[0]),
            max(a[1], b[1]),
            min(a[2], b[2]),
            max(a[3], b[3]),
            min(a[4], b[4]),
            max(a[5], b[5]),
        )
    # intersection
    return (
        max(a[0], b[0]),
        min(a[1], b[1]),
        max(a[2], b[2]),
        min(a[3], b[3]),
        max(a[4], b[4]),
        min(a[5], b[5]),
    )


class ImplicitShape(Shape):
    """Implicit shape defined by a signed distance / scalar field.

    The convention is ``f(x, y, z) < 0`` means inside the shape.

    :param func: scalar field ``(x, y, z) -> array``, or ``None`` for subclasses
    :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)`` or ``None``
    :param center: center of the shape
    :param orientation: orientation of the shape
    """

    def __init__(
        self: ImplicitShape,
        func: Field | None = None,
        bounds: BoundsType | None = None,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the implicit shape."""
        super().__init__(**kwargs)
        self._func = func
        self._bounds = bounds

    def evaluate(
        self: ImplicitShape,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Evaluate the scalar field at the given coordinates.

        :param x: x coordinates
        :param y: y coordinates
        :param z: z coordinates
        :return: scalar field values (negative = inside)
        """
        if self._func is None:
            err_msg = "No scalar field function defined"
            raise ValueError(err_msg)
        return self._func(x, y, z)

    def generate_vtk(
        self: ImplicitShape,
        bounds: BoundsType | None = None,
        resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate a VTK mesh of the implicit shape using the f < 0 convention.

        :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)``
        :param resolution: number of grid points per axis
        :return: triangulated surface mesh
        """
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

        clipped = grid.clip_scalar(scalars="implicit", invert=True)
        if clipped.n_cells == 0:
            return pv.PolyData()

        polydata = clipped.extract_surface().clean().triangulate()

        polydata = rotate_mesh(polydata, center=(0, 0, 0), rotation=self.orientation)
        return polydata.translate(xyz=self.center)

    def generate(
        self: ImplicitShape,
        bounds: BoundsType | None = None,
        resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> "cq.Shape":  # noqa: F821
        """Generate a CadQuery Shape from the implicit surface.

        :param bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)``
        :param resolution: number of grid points per axis
        :return: CadQuery Shape
        """
        import cadquery as cq

        from .tpms import ShellCreationError

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
                for start, end in zip(tri[:], tri[1:])
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

    def generateVtk(  # noqa: N802
        self: ImplicitShape,
        **kwargs: KwargsGenerateType,
    ) -> pv.PolyData:
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(**kwargs)

    # --- Boolean operators ---

    def __or__(self: ImplicitShape, other: ImplicitShape) -> ImplicitShape:
        """Union (``a | b``): inside where either field is negative."""
        return union(self, other)

    def __and__(self: ImplicitShape, other: ImplicitShape) -> ImplicitShape:
        """Intersection (``a & b``): inside where both fields are negative."""
        return intersection(self, other)

    def __sub__(self: ImplicitShape, other: ImplicitShape) -> ImplicitShape:
        """Difference (``a - b``): inside a but not b."""
        return difference(self, other)

    def __invert__(self: ImplicitShape) -> ImplicitShape:
        """Complement (``~a``): negate the field."""
        f = self._func
        return ImplicitShape(
            func=lambda x, y, z, _f=f: -_f(x, y, z),
            bounds=self._bounds,
        )

    # --- Smooth booleans ---

    def smooth_union(
        self: ImplicitShape,
        other: ImplicitShape,
        k: float,
    ) -> ImplicitShape:
        """Smooth union with blending radius *k*."""
        return smooth_union(self, other, k)

    def smooth_intersection(
        self: ImplicitShape,
        other: ImplicitShape,
        k: float,
    ) -> ImplicitShape:
        """Smooth intersection with blending radius *k*."""
        return smooth_intersection(self, other, k)

    def smooth_difference(
        self: ImplicitShape,
        other: ImplicitShape,
        k: float,
    ) -> ImplicitShape:
        """Smooth difference with blending radius *k*."""
        return smooth_difference(self, other, k)

    # --- Transforms ---

    def translate(  # type: ignore[override]
        self: ImplicitShape,
        offset: tuple[float, float, float],
    ) -> ImplicitShape:
        """Return a new shape translated by *offset*."""
        dx, dy, dz = offset
        f = self._func
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
        return ImplicitShape(
            func=lambda x, y, z, _f=f, _dx=dx, _dy=dy, _dz=dz: _f(
                x - _dx, y - _dy, z - _dz
            ),
            bounds=new_bounds,
        )

    def rotate(  # type: ignore[override]
        self: ImplicitShape,
        angles: tuple[float, float, float],
        convention: str = "ZXZ",
    ) -> ImplicitShape:
        """Return a new shape rotated by Euler *angles* (degrees)."""
        rot = Rotation.from_euler(convention, angles, degrees=True)
        inv_matrix = rot.inv().as_matrix()
        f = self._func
        return ImplicitShape(
            func=lambda x, y, z, _f=f, _m=inv_matrix: _f(
                *(_m @ np.array([x, y, z]))
            ),
            bounds=self._bounds,  # conservative: keep original bounds
        )

    def scale(
        self: ImplicitShape,
        factor: float,
    ) -> ImplicitShape:
        """Return a new shape uniformly scaled by *factor*."""
        f = self._func
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
        return ImplicitShape(
            func=lambda x, y, z, _f=f, _s=factor: _f(x / _s, y / _s, z / _s) * _s,
            bounds=new_bounds,
        )


# ---------------------------------------------------------------------------
# Module-level boolean operations
# ---------------------------------------------------------------------------


def union(a: ImplicitShape, b: ImplicitShape) -> ImplicitShape:
    """Union of two implicit shapes (hard boolean)."""
    fa, fb = a._func, b._func  # noqa: SLF001
    return ImplicitShape(
        func=lambda x, y, z, _fa=fa, _fb=fb: np.minimum(_fa(x, y, z), _fb(x, y, z)),
        bounds=_merge_bounds(a._bounds, b._bounds, "union"),  # noqa: SLF001
    )


def intersection(a: ImplicitShape, b: ImplicitShape) -> ImplicitShape:
    """Intersection of two implicit shapes (hard boolean)."""
    fa, fb = a._func, b._func  # noqa: SLF001
    return ImplicitShape(
        func=lambda x, y, z, _fa=fa, _fb=fb: np.maximum(_fa(x, y, z), _fb(x, y, z)),
        bounds=_merge_bounds(a._bounds, b._bounds, "intersection"),  # noqa: SLF001
    )


def difference(a: ImplicitShape, b: ImplicitShape) -> ImplicitShape:
    """Difference of two implicit shapes (a minus b)."""
    fa, fb = a._func, b._func  # noqa: SLF001
    return ImplicitShape(
        func=lambda x, y, z, _fa=fa, _fb=fb: np.maximum(
            _fa(x, y, z), -_fb(x, y, z)
        ),
        bounds=a._bounds,  # noqa: SLF001
    )


def smooth_union(
    a: ImplicitShape,
    b: ImplicitShape,
    k: float,
) -> ImplicitShape:
    """Smooth union with blending radius *k*."""
    fa, fb = a._func, b._func  # noqa: SLF001
    return ImplicitShape(
        func=lambda x, y, z, _fa=fa, _fb=fb, _k=k: _smooth_min(
            _fa(x, y, z), _fb(x, y, z), _k
        ),
        bounds=_merge_bounds(a._bounds, b._bounds, "union"),  # noqa: SLF001
    )


def smooth_intersection(
    a: ImplicitShape,
    b: ImplicitShape,
    k: float,
) -> ImplicitShape:
    """Smooth intersection with blending radius *k*."""
    fa, fb = a._func, b._func  # noqa: SLF001
    return ImplicitShape(
        func=lambda x, y, z, _fa=fa, _fb=fb, _k=k: _smooth_max(
            _fa(x, y, z), _fb(x, y, z), _k
        ),
        bounds=_merge_bounds(a._bounds, b._bounds, "intersection"),  # noqa: SLF001
    )


def smooth_difference(
    a: ImplicitShape,
    b: ImplicitShape,
    k: float,
) -> ImplicitShape:
    """Smooth difference (a minus b) with blending radius *k*."""
    fa, fb = a._func, b._func  # noqa: SLF001
    return ImplicitShape(
        func=lambda x, y, z, _fa=fa, _fb=fb, _k=k: _smooth_max(
            _fa(x, y, z), -_fb(x, y, z), _k
        ),
        bounds=a._bounds,  # noqa: SLF001
    )


def batch_smooth_union(
    shapes: list[ImplicitShape],
    k: float = 0.0,
) -> ImplicitShape:
    """Combine many shapes with smooth union in a flat loop (no recursion).

    This avoids the recursion-depth limit that arises when chaining hundreds
    of binary ``smooth_union`` calls, each wrapping the previous in a lambda.
    """
    if not shapes:
        msg = "batch_smooth_union requires at least one shape"
        raise ValueError(msg)

    funcs = [s._func for s in shapes]  # noqa: SLF001

    # Merge all bounds
    merged = shapes[0]._bounds  # noqa: SLF001
    for s in shapes[1:]:
        merged = _merge_bounds(merged, s._bounds, "union")  # noqa: SLF001

    def _batched(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        result = funcs[0](x, y, z)
        if k > 0:
            for fn in funcs[1:]:
                result = _smooth_min(result, fn(x, y, z), k)
        else:
            for fn in funcs[1:]:
                result = np.minimum(result, fn(x, y, z))
        return result

    return ImplicitShape(func=_batched, bounds=merged)


# ---------------------------------------------------------------------------
# Primitive factory functions
# ---------------------------------------------------------------------------


def implicit_sphere(
    center: tuple[float, float, float] = (0, 0, 0),
    radius: float = 1.0,
) -> ImplicitShape:
    """Sphere SDF: ``|p - c| - r``."""
    cx, cy, cz = center
    r = radius
    margin = r * 1.1
    return ImplicitShape(
        func=lambda x, y, z, _cx=cx, _cy=cy, _cz=cz, _r=r: np.sqrt(
            (x - _cx) ** 2 + (y - _cy) ** 2 + (z - _cz) ** 2
        )
        - _r,
        bounds=(
            cx - margin,
            cx + margin,
            cy - margin,
            cy + margin,
            cz - margin,
            cz + margin,
        ),
    )


def implicit_box(
    center: tuple[float, float, float] = (0, 0, 0),
    half_extents: tuple[float, float, float] = (0.5, 0.5, 0.5),
) -> ImplicitShape:
    """Axis-aligned box SDF (exact)."""
    cx, cy, cz = center
    hx, hy, hz = half_extents
    margin = max(hx, hy, hz) * 0.1

    def _box_sdf(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        qx = np.abs(x - cx) - hx
        qy = np.abs(y - cy) - hy
        qz = np.abs(z - cz) - hz
        outside = np.sqrt(
            np.maximum(qx, 0.0) ** 2
            + np.maximum(qy, 0.0) ** 2
            + np.maximum(qz, 0.0) ** 2
        )
        inside = np.minimum(np.maximum(qx, np.maximum(qy, qz)), 0.0)
        return outside + inside

    return ImplicitShape(
        func=_box_sdf,
        bounds=(
            cx - hx - margin,
            cx + hx + margin,
            cy - hy - margin,
            cy + hy + margin,
            cz - hz - margin,
            cz + hz + margin,
        ),
    )


def implicit_cylinder(
    center: tuple[float, float, float] = (0, 0, 0),
    axis: tuple[float, float, float] = (0, 0, 1),
    radius: float = 0.5,
    height: float = 1.0,
) -> ImplicitShape:
    """Capped cylinder SDF along arbitrary axis."""
    c = np.asarray(center, dtype=np.float64)
    a = np.asarray(axis, dtype=np.float64)
    a = a / np.linalg.norm(a)
    half_h = height / 2.0

    def _cyl_sdf(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        px, py, pz = x - c[0], y - c[1], z - c[2]
        proj = px * a[0] + py * a[1] + pz * a[2]
        perp_x = px - proj * a[0]
        perp_y = py - proj * a[1]
        perp_z = pz - proj * a[2]
        d_radial = np.sqrt(perp_x**2 + perp_y**2 + perp_z**2) - radius
        d_axial = np.abs(proj) - half_h
        outside = np.sqrt(
            np.maximum(d_radial, 0.0) ** 2 + np.maximum(d_axial, 0.0) ** 2
        )
        inside = np.minimum(np.maximum(d_radial, d_axial), 0.0)
        return outside + inside

    extent = max(radius, half_h) * 1.1 + np.linalg.norm(c)
    return ImplicitShape(
        func=_cyl_sdf,
        bounds=(
            c[0] - extent,
            c[0] + extent,
            c[1] - extent,
            c[1] + extent,
            c[2] - extent,
            c[2] + extent,
        ),
    )


def implicit_plane(
    point: tuple[float, float, float] = (0, 0, 0),
    normal: tuple[float, float, float] = (0, 0, 1),
) -> ImplicitShape:
    """Half-space SDF: ``dot(p - point, normal)``."""
    pt = np.asarray(point, dtype=np.float64)
    n = np.asarray(normal, dtype=np.float64)
    n = n / np.linalg.norm(n)

    return ImplicitShape(
        func=lambda x, y, z: (x - pt[0]) * n[0]
        + (y - pt[1]) * n[1]
        + (z - pt[2]) * n[2],
        bounds=None,
    )


def implicit_torus(
    center: tuple[float, float, float] = (0, 0, 0),
    major_r: float = 1.0,
    minor_r: float = 0.25,
) -> ImplicitShape:
    """Torus SDF in XZ plane."""
    cx, cy, cz = center

    def _torus_sdf(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        px, py, pz = x - cx, y - cy, z - cz
        q_xz = np.sqrt(px**2 + pz**2) - major_r
        return np.sqrt(q_xz**2 + py**2) - minor_r

    margin = (major_r + minor_r) * 1.1
    return ImplicitShape(
        func=_torus_sdf,
        bounds=(
            cx - margin,
            cx + margin,
            cy - margin,
            cy + margin,
            cz - margin,
            cz + margin,
        ),
    )


def implicit_capsule(
    start: tuple[float, float, float] = (0, 0, -0.5),
    end: tuple[float, float, float] = (0, 0, 0.5),
    radius: float = 0.25,
) -> ImplicitShape:
    """Capsule SDF (line segment distance minus radius)."""
    a = np.asarray(start, dtype=np.float64)
    b = np.asarray(end, dtype=np.float64)
    ab = b - a
    ab_dot = float(np.dot(ab, ab))

    if ab_dot < 1e-30:
        # Degenerate capsule (start == end): fall back to sphere
        return implicit_sphere(center=tuple(a), radius=radius)

    def _capsule_sdf(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        px, py, pz = x - a[0], y - a[1], z - a[2]
        t = np.clip((px * ab[0] + py * ab[1] + pz * ab[2]) / ab_dot, 0.0, 1.0)
        dx = px - t * ab[0]
        dy = py - t * ab[1]
        dz = pz - t * ab[2]
        return np.sqrt(dx**2 + dy**2 + dz**2) - radius

    mid = (a + b) / 2.0
    extent = np.linalg.norm(ab) / 2.0 + radius
    margin = extent * 1.1
    return ImplicitShape(
        func=_capsule_sdf,
        bounds=(
            mid[0] - margin,
            mid[0] + margin,
            mid[1] - margin,
            mid[1] + margin,
            mid[2] - margin,
            mid[2] + margin,
        ),
    )


def implicit_ellipsoid(
    center: tuple[float, float, float] = (0, 0, 0),
    radii: tuple[float, float, float] = (1.0, 0.5, 0.25),
) -> ImplicitShape:
    """Ellipsoid implicit field: ``sqrt((x/rx)^2 + (y/ry)^2 + (z/rz)^2) - 1``.

    This is an approximate SDF but has the correct zero-level set.
    """
    cx, cy, cz = center
    rx, ry, rz = radii
    margin = max(rx, ry, rz) * 1.1

    def _ellipsoid_field(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return (
            np.sqrt(
                ((x - cx) / rx) ** 2
                + ((y - cy) / ry) ** 2
                + ((z - cz) / rz) ** 2
            )
            - 1.0
        )

    return ImplicitShape(
        func=_ellipsoid_field,
        bounds=(
            cx - margin,
            cx + margin,
            cy - margin,
            cy + margin,
            cz - margin,
            cz + margin,
        ),
    )


# ---------------------------------------------------------------------------
# Utility operations
# ---------------------------------------------------------------------------


def shell(shape: ImplicitShape, thickness: float) -> ImplicitShape:
    """Hollow shell: ``|f(p)| - thickness / 2``."""
    f = shape._func  # noqa: SLF001
    half_t = thickness / 2.0
    return ImplicitShape(
        func=lambda x, y, z, _f=f, _ht=half_t: np.abs(_f(x, y, z)) - _ht,
        bounds=shape._bounds,  # noqa: SLF001
    )


def repeat(
    shape: ImplicitShape,
    spacing: tuple[float, float, float],
    k: float = 0.0,
) -> ImplicitShape:
    """Infinite repetition via coordinate modulo.

    :param shape: unit cell shape to tile
    :param spacing: ``(sx, sy, sz)`` repetition period per axis
    :param k: smooth blending radius across cell boundaries.
        When ``k > 0``, the base field is evaluated at the 26 neighboring
        periodic images in addition to the current cell and all values
        are combined with smooth minimum, so that primitives from adjacent
        cells blend seamlessly.  When ``k <= 0`` (default), a simple
        coordinate-modulo repetition is used (hard tiling).
    """
    sx, sy, sz = spacing
    f = shape._func  # noqa: SLF001

    if k <= 0:

        def _repeated(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            rx = np.mod(x + sx / 2, sx) - sx / 2
            ry = np.mod(y + sy / 2, sy) - sy / 2
            rz = np.mod(z + sz / 2, sz) - sz / 2
            return f(rx, ry, rz)

    else:
        offsets = [
            (dx * sx, dy * sy, dz * sz)
            for dx in (-1, 0, 1)
            for dy in (-1, 0, 1)
            for dz in (-1, 0, 1)
        ]

        def _repeated(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            # Wrap to the central cell
            cx = np.mod(x + sx / 2, sx) - sx / 2
            cy = np.mod(y + sy / 2, sy) - sy / 2
            cz = np.mod(z + sz / 2, sz) - sz / 2
            # Evaluate at all 27 images, accumulate with smooth min
            result = f(cx + offsets[0][0], cy + offsets[0][1], cz + offsets[0][2])
            for ox, oy, oz in offsets[1:]:
                result = _smooth_min(result, f(cx + ox, cy + oy, cz + oz), k)
            return result

    return ImplicitShape(func=_repeated, bounds=None)


def blend(
    a: ImplicitShape,
    b: ImplicitShape,
    factor: float = 0.5,
) -> ImplicitShape:
    """Linear interpolation between two fields: ``(1-t)*a + t*b``."""
    fa, fb = a._func, b._func  # noqa: SLF001
    t = factor
    return ImplicitShape(
        func=lambda x, y, z, _fa=fa, _fb=fb, _t=t: (1.0 - _t) * _fa(x, y, z)
        + _t * _fb(x, y, z),
        bounds=_merge_bounds(a._bounds, b._bounds, "union"),  # noqa: SLF001
    )


def from_field(func: Field, bounds: BoundsType | None = None) -> ImplicitShape:
    """Wrap any callable ``f(x, y, z) -> scalar`` as an ImplicitShape."""
    return ImplicitShape(func=func, bounds=bounds)
