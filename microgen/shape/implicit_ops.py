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

from . import shape as _shape

if TYPE_CHECKING:
    from .shape import BoundsType, Field, Shape


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _make_shape(func: Field, bounds: BoundsType | None) -> Shape:
    """Create a bare Shape with only an implicit scalar field."""
    return _shape.Shape(func=func, bounds=bounds)


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
            _fa(x, y, z),
            -_fb(x, y, z),
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
            _fa(x, y, z),
            _fb(x, y, z),
            _k,
        ),
        bounds=_merge_bounds(a.bounds, b.bounds, "union"),
    )


def smooth_intersection(a: Shape, b: Shape, k: float) -> Shape:
    """Smooth intersection with blending radius *k*."""
    fa, fb = a.require_func(), b.require_func()
    return _make_shape(
        func=lambda x, y, z, _fa=fa, _fb=fb, _k=k: _smooth_max(
            _fa(x, y, z),
            _fb(x, y, z),
            _k,
        ),
        bounds=_merge_bounds(a.bounds, b.bounds, "intersection"),
    )


def smooth_difference(a: Shape, b: Shape, k: float) -> Shape:
    """Smooth difference (a minus b) with blending radius *k*."""
    fa, fb = a.require_func(), b.require_func()
    return _make_shape(
        func=lambda x, y, z, _fa=fa, _fb=fb, _k=k: _smooth_max(
            _fa(x, y, z),
            -_fb(x, y, z),
            _k,
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


def shell(shape: Shape, thickness: float | Field) -> Shape:
    """Hollow shell: ``|f(p)| - thickness(p) / 2``.

    ``thickness`` may be a constant (scalar) or a callable
    ``thickness(x, y, z) -> array`` for spatially-varying shells.  Negative or
    zero thickness at a point yields no inclusion in the shell at that point.
    """
    f = shape.require_func()
    if callable(thickness):
        t_func = thickness

        def _shell_field(x, y, z, _f=f, _t=t_func):  # noqa: ANN001
            return np.abs(_f(x, y, z)) - _t(x, y, z) / 2.0

        return _make_shape(func=_shell_field, bounds=shape.bounds)

    half_t = float(thickness) / 2.0
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


def box(
    dims: tuple[float, float, float],
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> Shape:
    """Axis-aligned box as an F-rep Shape.

    SDF formula ``max(|x-cx|-hx, |y-cy|-hy, |z-cz|-hz)``: signed distance to
    the box surface (negative inside, positive outside, zero on the surface).
    Useful as a clipping primitive — e.g. ``intersection(skeletal, box(...))``
    bounds an unbounded TPMS skeletal field to a single cell.

    :param dims: full side lengths ``(dx, dy, dz)``
    :param center: box center (default origin)
    :return: :class:`~microgen.shape.shape.Shape` carrying the box SDF
    """
    hx, hy, hz = (0.5 * float(d) for d in dims)
    cx, cy, cz = (float(c) for c in center)

    def _box_sdf(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return np.maximum.reduce(
            [np.abs(x - cx) - hx, np.abs(y - cy) - hy, np.abs(z - cz) - hz],
        )

    bounds: BoundsType = (cx - hx, cx + hx, cy - hy, cy + hy, cz - hz, cz + hz)
    return _make_shape(func=_box_sdf, bounds=bounds)


def _fd_sdf(
    f: Field,
    epsilon: float,
) -> Field:
    """SDF via central finite differences (fallback)."""

    def sdf(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        h = 1e-5
        val = f(x, y, z)
        gx = (f(x + h, y, z) - f(x - h, y, z)) / (2 * h)
        gy = (f(x, y + h, z) - f(x, y - h, z)) / (2 * h)
        gz = (f(x, y, z + h) - f(x, y, z - h)) / (2 * h)
        grad_mag = np.sqrt(gx**2 + gy**2 + gz**2)
        # Where the gradient is degenerate (e.g. flat-z fields like the
        # honeycomb_* surfaces, or saddle points), normalization would
        # blow up to ±1/epsilon — preserve the raw field's *sign* by
        # falling back to the unnormalized value there.
        return np.where(grad_mag > epsilon, val / np.maximum(grad_mag, epsilon), val)

    return sdf


def normalize_to_sdf(shape: Shape, epsilon: float = 1e-10) -> Shape:
    """Return a new Shape with gradient-normalized SDF field: ``f / |nabla f|``.

    Uses ``autograd`` for exact analytical gradients when the field function
    is differentiable through ``autograd.numpy``.  Falls back to central
    finite differences otherwise.

    :param shape: shape whose implicit field to normalize
    :param epsilon: floor for gradient magnitude (avoids division by zero
        at saddle points)
    """
    f = shape.require_func()

    # Try autograd first; fall back to FD if it fails at construction
    # OR at first evaluation (autograd may succeed at construction but
    # fail when the inner function uses non-autograd numpy ops).
    try:
        from autograd import elementwise_grad  # noqa: PLC0415

        dfdx = elementwise_grad(f, argnum=0)
        dfdy = elementwise_grad(f, argnum=1)
        dfdz = elementwise_grad(f, argnum=2)

        # Probe all three gradients to detect deferred failures
        # (autograd wraps only the argnum-th arg, so a non-autograd op
        # on y/z would pass the dfdx probe but crash in dfdy/dfdz)
        _probe = np.array([0.0])
        dfdx(_probe, _probe, _probe)
        dfdy(_probe, _probe, _probe)
        dfdz(_probe, _probe, _probe)

        def sdf(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            val = f(x, y, z)
            grad_mag = np.sqrt(
                dfdx(x, y, z) ** 2 + dfdy(x, y, z) ** 2 + dfdz(x, y, z) ** 2,
            )
            # Same fallback as in `_fd_sdf`: where the gradient vanishes
            # (degenerate flat-z fields, saddle points), keep the raw value
            # so its sign is preserved without exploding into ±1/epsilon.
            return np.where(
                grad_mag > epsilon, val / np.maximum(grad_mag, epsilon), val
            )

    except Exception:  # noqa: BLE001
        sdf = _fd_sdf(f, epsilon)

    return _make_shape(func=sdf, bounds=shape.bounds)


def variable_shell(
    shape: Shape,
    thickness_func: Field,
) -> Shape:
    """Shell with spatially-varying thickness: ``|f(p)| - t(p)/2``.

    :param shape: shape whose implicit field defines the surface
    :param thickness_func: callable ``(x, y, z) -> thickness`` returning
        the local shell thickness
    """
    f = shape.require_func()

    def _var_shell(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return np.abs(f(x, y, z)) - thickness_func(x, y, z) / 2.0

    return _make_shape(func=_var_shell, bounds=shape.bounds)
