"""F-rep Implicit Operations.

==========================================================
Implicit Operations (:mod:`microgen.shape.implicit_ops`)
==========================================================

Module-level boolean, blending, and utility operations for shapes
that carry an implicit scalar field (``_func``).  All functions accept
and return :class:`~microgen.shape.shape.Shape` instances.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from .shape import BoundsType, Field, Shape


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _make_shape(func: Field, bounds: BoundsType | None) -> Shape:
    """Create a Shape with an implicit field (single deferred import)."""
    from .shape import Shape  # noqa: PLC0415

    return Shape(func=func, bounds=bounds)


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


# ---------------------------------------------------------------------------
# Unary operations
# ---------------------------------------------------------------------------


def complement(a: Shape) -> Shape:
    """Complement (negate the field): inside becomes outside and vice versa."""
    f = a.require_func()
    return _make_shape(
        func=lambda x, y, z, _f=f: -_f(x, y, z),
        bounds=a.bounds,
    )


# ---------------------------------------------------------------------------
# Hard boolean operations
# ---------------------------------------------------------------------------


def union(a: Shape, b: Shape) -> Shape:
    """Union of two shapes (hard boolean)."""
    fa, fb = a.require_func(), b.require_func()
    return _make_shape(
        func=lambda x, y, z, _fa=fa, _fb=fb: np.minimum(_fa(x, y, z), _fb(x, y, z)),
        bounds=_merge_bounds(a.bounds, b.bounds, "union"),
    )


def intersection(a: Shape, b: Shape) -> Shape:
    """Intersection of two shapes (hard boolean)."""
    fa, fb = a.require_func(), b.require_func()
    return _make_shape(
        func=lambda x, y, z, _fa=fa, _fb=fb: np.maximum(_fa(x, y, z), _fb(x, y, z)),
        bounds=_merge_bounds(a.bounds, b.bounds, "intersection"),
    )


def difference(a: Shape, b: Shape) -> Shape:
    """Difference of two shapes (a minus b)."""
    fa, fb = a.require_func(), b.require_func()
    return _make_shape(
        func=lambda x, y, z, _fa=fa, _fb=fb: np.maximum(
            _fa(x, y, z), -_fb(x, y, z),
        ),
        bounds=a.bounds,
    )


# ---------------------------------------------------------------------------
# Smooth boolean operations
# ---------------------------------------------------------------------------


def smooth_union(a: Shape, b: Shape, k: float) -> Shape:
    """Smooth union with blending radius *k*."""
    fa, fb = a.require_func(), b.require_func()
    return _make_shape(
        func=lambda x, y, z, _fa=fa, _fb=fb, _k=k: _smooth_min(
            _fa(x, y, z), _fb(x, y, z), _k,
        ),
        bounds=_merge_bounds(a.bounds, b.bounds, "union"),
    )


def smooth_intersection(a: Shape, b: Shape, k: float) -> Shape:
    """Smooth intersection with blending radius *k*."""
    fa, fb = a.require_func(), b.require_func()
    return _make_shape(
        func=lambda x, y, z, _fa=fa, _fb=fb, _k=k: _smooth_max(
            _fa(x, y, z), _fb(x, y, z), _k,
        ),
        bounds=_merge_bounds(a.bounds, b.bounds, "intersection"),
    )


def smooth_difference(a: Shape, b: Shape, k: float) -> Shape:
    """Smooth difference (a minus b) with blending radius *k*."""
    fa, fb = a.require_func(), b.require_func()
    return _make_shape(
        func=lambda x, y, z, _fa=fa, _fb=fb, _k=k: _smooth_max(
            _fa(x, y, z), -_fb(x, y, z), _k,
        ),
        bounds=a.bounds,
    )


# ---------------------------------------------------------------------------
# Batch operations
# ---------------------------------------------------------------------------


def batch_smooth_union(
    shapes: list[Shape],
    k: float = 0.0,
) -> Shape:
    """Combine many shapes with smooth union in a flat loop (no recursion).

    This avoids the recursion-depth limit that arises when chaining hundreds
    of binary ``smooth_union`` calls, each wrapping the previous in a lambda.
    """
    if not shapes:
        msg = "batch_smooth_union requires at least one shape"
        raise ValueError(msg)

    funcs = [s.require_func() for s in shapes]

    merged = shapes[0].bounds
    for s in shapes[1:]:
        merged = _merge_bounds(merged, s.bounds, "union")

    def _batched(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        if k > 0:
            result = funcs[0](x, y, z)
            for fn in funcs[1:]:
                result = _smooth_min(result, fn(x, y, z), k)
            return result
        # Hard union: vectorized reduction
        all_fields = np.stack([fn(x, y, z) for fn in funcs], axis=0)
        return np.min(all_fields, axis=0)

    return _make_shape(func=_batched, bounds=merged)


# ---------------------------------------------------------------------------
# Utility operations
# ---------------------------------------------------------------------------


def shell(shape: Shape, thickness: float) -> Shape:
    """Hollow shell: ``|f(p)| - thickness / 2``."""
    f = shape.require_func()
    half_t = thickness / 2.0
    return _make_shape(
        func=lambda x, y, z, _f=f, _ht=half_t: np.abs(_f(x, y, z)) - _ht,
        bounds=shape.bounds,
    )


def repeat(
    shape: Shape,
    spacing: tuple[float, float, float],
    k: float = 0.0,
) -> Shape:
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
    f = shape.require_func()

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
            cx = np.mod(x + sx / 2, sx) - sx / 2
            cy = np.mod(y + sy / 2, sy) - sy / 2
            cz = np.mod(z + sz / 2, sz) - sz / 2
            result = f(cx + offsets[0][0], cy + offsets[0][1], cz + offsets[0][2])
            for ox, oy, oz in offsets[1:]:
                result = _smooth_min(result, f(cx + ox, cy + oy, cz + oz), k)
            return result

    return _make_shape(func=_repeated, bounds=None)


def blend(
    a: Shape,
    b: Shape,
    factor: float = 0.5,
) -> Shape:
    """Linear interpolation between two fields: ``(1-t)*a + t*b``."""
    fa, fb = a.require_func(), b.require_func()
    t = factor
    return _make_shape(
        func=lambda x, y, z, _fa=fa, _fb=fb, _t=t: (1.0 - _t) * _fa(x, y, z)
        + _t * _fb(x, y, z),
        bounds=_merge_bounds(a.bounds, b.bounds, "union"),
    )


def from_field(
    func: Field,
    bounds: BoundsType | None = None,
) -> Shape:
    """Wrap any callable ``f(x, y, z) -> scalar`` as a Shape with an implicit field."""
    return _make_shape(func=func, bounds=bounds)
