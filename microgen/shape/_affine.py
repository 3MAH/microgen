r"""Affine-transform helpers shared by :class:`Shape` and :class:`Phase`.

Every transform is the world map :math:`T(x) = A (x - p) + p + t` with
:math:`A` a rotation or a positive diagonal scale, :math:`p` the pivot and
:math:`t` a shift.  A field :math:`f` maps to
:math:`f_T(x) = c\, f(A^{-1}(x - p - t) + p)` with :math:`c > 0`; the solid
:math:`\{f < \mathrm{iso}\}` maps exactly to its image when ``iso`` is
multiplied by the same :math:`c`.  :math:`c = 1` for rigid maps and
:math:`c = \min_i s_i` for scales, so a signed distance stays exact under
similarities and becomes a 1-Lipschitz lower bound under per-axis scales.
The period lattice :math:`\mathrm{diag}(L)\,\mathbb{Z}^3` maps to
:math:`A\,\mathrm{diag}(L)\,\mathbb{Z}^3`.
"""

from __future__ import annotations

import itertools
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cached_property

import numpy as np
import numpy.typing as npt
from scipy.spatial.transform import Rotation

from ._types import BoundsType, Field, PeriodType

# Snapping tolerance of a rotation matrix to a signed permutation.
_SNAP_TOL = 1e-12

# Largest deviation of ``M^T M`` from the identity accepted for a rotation
# matrix; SciPy would otherwise replace any matrix by its nearest rotation.
_ORTHO_TOL = 1e-9


def signed_permutation(matrix: npt.NDArray[np.float64]) -> npt.NDArray | None:
    """Return the exact signed permutation ``matrix`` rounds to, or ``None``."""
    rounded = np.round(matrix)
    if np.max(np.abs(matrix - rounded)) > _SNAP_TOL:
        return None
    ones = np.ones(3)
    if not (
        np.array_equal(np.abs(rounded).sum(axis=0), ones)
        and np.array_equal(np.abs(rounded).sum(axis=1), ones)
    ):
        return None
    return rounded


def rotation_matrix(rotation: Rotation) -> npt.NDArray[np.float64]:
    """Return the matrix of ``rotation``, exact when it is a signed permutation.

    SciPy stores rotations as quaternions, so ``as_matrix`` of a quarter
    turn carries ``1 +- 2e-16`` entries; snapping restores the exact matrix.
    """
    matrix = rotation.as_matrix()
    snapped = signed_permutation(matrix)
    return snapped if snapped is not None else matrix


def as_rotation(rotation: Rotation | npt.ArrayLike) -> Rotation:
    """Coerce to a SciPy ``Rotation``.

    :param rotation: a ``Rotation`` or a proper orthogonal 3x3 matrix
    :return: the rotation
    :raises ValueError: if the matrix is not orthogonal or has ``det < 0``
    """
    if isinstance(rotation, Rotation):
        return rotation
    matrix = np.asarray(rotation, dtype=np.float64)
    if matrix.shape != (3, 3) or (
        np.max(np.abs(matrix.T @ matrix - np.eye(3))) > _ORTHO_TOL
    ):
        err_msg = f"rotation matrix must be orthogonal. Given: {matrix.tolist()}"
        raise ValueError(err_msg)
    return Rotation.from_matrix(matrix)


def scale_factors(factor: float | Sequence[float]) -> npt.NDArray[np.float64]:
    """Return ``factor`` as a positive ``(sx, sy, sz)`` array.

    :raises ValueError: if any component is not strictly positive
        (reflections are not supported).
    """
    s = np.broadcast_to(np.asarray(factor, dtype=np.float64), (3,)).copy()
    if np.any(s <= 0.0):
        err_msg = f"scale factors must be strictly positive. Given: {factor}"
        raise ValueError(err_msg)
    return s


def transform_bounds(
    bounds: BoundsType,
    matrix: npt.NDArray[np.float64],
    point: npt.NDArray[np.float64],
    shift: npt.NDArray[np.float64],
) -> BoundsType:
    """Return the AABB of the image of the 8 corners of ``bounds``."""
    corners = np.array(list(itertools.product(bounds[0:2], bounds[2:4], bounds[4:6])))
    mapped = (corners - point) @ matrix.T + point + shift
    lo = mapped.min(axis=0)
    hi = mapped.max(axis=0)
    return (
        float(lo[0]),
        float(hi[0]),
        float(lo[1]),
        float(hi[1]),
        float(lo[2]),
        float(hi[2]),
    )


def transform_field(
    field: Field,
    matrix: npt.NDArray[np.float64],
    point: npt.NDArray[np.float64],
    shift: npt.NDArray[np.float64],
    multiplier: float = 1.0,
) -> Field:
    """Return ``x -> multiplier * field(matrix^{-1} (x - point - shift) + point)``.

    Translations and diagonal scales get per-axis closures (3 operations
    per axis instead of a dense 3x3 product).
    """
    px, py, pz = (float(v) for v in point)
    qx, qy, qz = (float(v) for v in point + shift)
    c = float(multiplier)

    if np.count_nonzero(matrix - np.diag(np.diag(matrix))) == 0:
        ix, iy, iz = (1.0 / float(v) for v in np.diag(matrix))

        def _diagonal(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            values = field((x - qx) * ix + px, (y - qy) * iy + py, (z - qz) * iz + pz)
            return values if c == 1.0 else c * values

        return _diagonal

    m = np.linalg.inv(matrix)

    def _field(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        xx = x - qx
        yy = y - qy
        zz = z - qz
        lx = m[0, 0] * xx + m[0, 1] * yy + m[0, 2] * zz + px
        ly = m[1, 0] * xx + m[1, 1] * yy + m[1, 2] * zz + py
        lz = m[2, 0] * xx + m[2, 1] * yy + m[2, 2] * zz + pz
        values = field(lx, ly, lz)
        return values if c == 1.0 else c * values

    return _field


def rotate_period(
    period: PeriodType | None,
    matrix: npt.NDArray[np.float64],
) -> PeriodType | None:
    """Period after a rotation: permuted for a signed permutation, else ``None``.

    The lattice ``diag(L) Z^3`` maps to ``R diag(L) Z^3``, which is
    axis-aligned only when ``R`` is a signed permutation.
    """
    if period is None:
        return None
    perm = signed_permutation(matrix)
    if perm is None:
        return None
    rotated = np.abs(perm) @ np.asarray(period, dtype=np.float64)
    return (float(rotated[0]), float(rotated[1]), float(rotated[2]))


@dataclass(frozen=True)
class AffineMap:
    """The map ``x -> matrix (x - point) + point + shift`` and its field factor.

    Built by :func:`translation`, :func:`rotation_about` or :func:`scaling`;
    ``rotation`` and ``factors`` record which one, for the period and the
    native-parameter updates.
    """

    matrix: npt.NDArray[np.float64]
    point: npt.NDArray[np.float64]
    shift: npt.NDArray[np.float64]
    multiplier: float = 1.0
    rotation: Rotation | None = None
    factors: npt.NDArray[np.float64] | None = None

    @cached_property
    def offset(self) -> npt.NDArray[np.float64]:
        """Translation part ``b`` of the same map written ``x -> matrix x + b``."""
        return self.point + self.shift - self.matrix @ self.point

    @cached_property
    def homogeneous(self) -> npt.NDArray[np.float64]:
        """The ``(4, 4)`` homogeneous matrix (PyVista ``transform``)."""
        out = np.eye(4)
        out[:3, :3] = self.matrix
        out[:3, 3] = self.offset
        return out

    def apply_point(self, xyz: Sequence[float]) -> tuple[float, float, float]:
        """Map one point."""
        out = self.matrix @ np.asarray(xyz, dtype=np.float64) + self.offset
        return (float(out[0]), float(out[1]), float(out[2]))

    def apply_bounds(self, bounds: BoundsType) -> BoundsType:
        """AABB of the mapped ``bounds``."""
        return transform_bounds(bounds, self.matrix, self.point, self.shift)

    def apply_field(self, field: Field) -> Field:
        """The mapped field, multiplied by :attr:`multiplier`."""
        return transform_field(
            field, self.matrix, self.point, self.shift, self.multiplier
        )

    def apply_period(self, period: PeriodType | None) -> PeriodType | None:
        """Period after the map: kept, permuted or reset, or scaled."""
        if period is None:
            return None
        if self.rotation is not None:
            return rotate_period(period, self.matrix)
        if self.factors is not None:
            return (
                float(self.factors[0] * period[0]),
                float(self.factors[1] * period[1]),
                float(self.factors[2] * period[2]),
            )
        return period

    @cached_property
    def grows_bounds(self) -> bool:
        """``True`` for a rotation whose image AABB is larger than the image box."""
        return self.rotation is not None and signed_permutation(self.matrix) is None


def _as_point(point: Sequence[float] | None) -> npt.NDArray[np.float64]:
    if point is None:
        return np.zeros(3)
    return np.asarray(point, dtype=np.float64).reshape(3)


def translation(offset: Sequence[float]) -> AffineMap:
    """Translation by ``offset``."""
    shift = np.asarray(offset, dtype=np.float64).reshape(3)
    return AffineMap(np.eye(3), np.zeros(3), shift)


def rotation_about(
    rotation: Rotation | npt.ArrayLike,
    point: Sequence[float] | None = None,
) -> AffineMap:
    """Rotation about ``point`` (world origin when ``None``)."""
    rot = as_rotation(rotation)
    return AffineMap(rotation_matrix(rot), _as_point(point), np.zeros(3), rotation=rot)


def scaling(
    factor: float | Sequence[float],
    point: Sequence[float] | None = None,
) -> AffineMap:
    """Positive diagonal scale about ``point``, with field factor ``min(s)``."""
    factors = scale_factors(factor)
    return AffineMap(
        np.diag(factors),
        _as_point(point),
        np.zeros(3),
        multiplier=float(factors.min()),
        factors=factors,
    )
