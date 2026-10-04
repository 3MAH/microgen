"""Affine-transform helpers shared by :class:`Shape` and :class:`Phase`.

Every transform is the world map ``x -> A (x - p) + p + t`` with ``A`` a
rotation or a positive diagonal scale, ``p`` the pivot and ``t`` a shift.
A field ``f`` maps to ``c * f(A^{-1} (x - p - t) + p)`` with ``c > 0``; the
solid ``{f < iso}`` maps exactly to its image when ``iso`` is multiplied by
the same ``c``.  ``c = 1`` for rigid maps and ``c = min(s)`` for scales, so
an SDF stays an SDF under similarities and a 1-Lipschitz lower bound under
anisotropic scales.
"""

from __future__ import annotations

import itertools
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt
from scipy.spatial.transform import Rotation

from ._types import BoundsType, Field, PeriodType

# Entries closer than this to 0 / ±1 are snapped when a rotation is a
# signed permutation (e.g. ``Rotation.from_euler("z", 90, degrees=True)``
# carries ``cos(pi/2) ~ 6e-17``), so periodic fields stay bit-exact.
_SNAP_TOL = 1e-12


def signed_permutation(matrix: npt.NDArray[np.float64]) -> npt.NDArray | None:
    """Return the exact signed permutation ``matrix`` rounds to, or ``None``."""
    rounded = np.round(matrix)
    if np.max(np.abs(matrix - rounded)) > _SNAP_TOL:
        return None
    if not np.array_equal(np.abs(rounded).sum(axis=0), np.ones(3)) or not (
        np.array_equal(np.abs(rounded).sum(axis=1), np.ones(3))
    ):
        return None
    return rounded


def as_rotation(rotation: Rotation | npt.ArrayLike) -> Rotation:
    """Coerce to a SciPy ``Rotation``, snapping near signed permutations.

    A 3x3 matrix goes through ``Rotation.from_matrix``, which rejects
    improper matrices (``det < 0``).
    """
    rot = (
        rotation
        if isinstance(rotation, Rotation)
        else Rotation.from_matrix(np.asarray(rotation, dtype=np.float64))
    )
    snapped = signed_permutation(rot.as_matrix())
    return Rotation.from_matrix(snapped) if snapped is not None else rot


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


def as_point(point: Sequence[float] | None) -> npt.NDArray[np.float64]:
    """Return the pivot as a float64 array, the world origin when ``None``."""
    if point is None:
        return np.zeros(3)
    return np.asarray(point, dtype=np.float64).reshape(3)


def transform_point(
    xyz: Sequence[float],
    matrix: npt.NDArray[np.float64],
    point: npt.NDArray[np.float64],
    shift: npt.NDArray[np.float64],
) -> tuple[float, float, float]:
    """Map a point through ``x -> matrix (x - point) + point + shift``."""
    out = matrix @ (np.asarray(xyz, dtype=np.float64) - point) + point + shift
    return (float(out[0]), float(out[1]), float(out[2]))


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
    """Return ``x -> multiplier * field(matrix^{-1} (x - point - shift) + point)``."""
    m = np.linalg.inv(matrix)
    px, py, pz = (float(v) for v in point)
    qx, qy, qz = (float(v) for v in point + shift)
    c = float(multiplier)

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
    rotation: Rotation,
) -> PeriodType | None:
    """Period after a rotation: permuted for a signed permutation, else ``None``.

    The period lattice ``diag(L) Z^3`` maps to ``R diag(L) Z^3``, which is
    axis-aligned only when ``R`` is a signed permutation.
    """
    if period is None:
        return None
    perm = signed_permutation(rotation.as_matrix())
    if perm is None:
        return None
    rotated = np.abs(perm) @ np.asarray(period, dtype=np.float64)
    return (float(rotated[0]), float(rotated[1]), float(rotated[2]))


def scale_period(
    period: PeriodType | None,
    factors: npt.NDArray[np.float64],
) -> PeriodType | None:
    """Period after a positive diagonal scale: ``(sx Lx, sy Ly, sz Lz)``."""
    if period is None:
        return None
    return (
        float(factors[0] * period[0]),
        float(factors[1] * period[1]),
        float(factors[2] * period[2]),
    )
