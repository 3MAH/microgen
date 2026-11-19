"""Deterministic interpolation of offsets sampled on a fixed TPMS grid."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import numpy.typing as npt

    from ._types import Field


def sampled_field(
    coordinates: npt.NDArray[np.float64], values: npt.NDArray[np.float64]
) -> Field:
    """Freeze samples; interpolate linearly, with a defined exterior extension.

    Tensor grids use trilinear interpolation and linear extrapolation. Clipped
    grids use simplex interpolation and nearest-sample extension outside their
    convex hull. The clipping envelope still determines the physical domain.
    """
    from scipy.interpolate import (  # noqa: PLC0415
        LinearNDInterpolator,
        NearestNDInterpolator,
        RegularGridInterpolator,
    )

    points = np.array(coordinates, dtype=np.float64, copy=True)
    samples = np.array(values, dtype=np.float64, copy=True).ravel(order="F")
    if samples.size != len(points) or not np.isfinite(samples).all():
        raise ValueError("offset samples must be finite and match the original grid")
    unique, first, inverse = np.unique(
        points, axis=0, return_index=True, return_inverse=True
    )
    if not np.allclose(samples, samples[first][inverse], rtol=1e-12, atol=1e-14):
        raise ValueError("offset samples conflict at coincident grid points")
    points, samples = unique, samples[first]
    axes = tuple(np.unique(points[:, axis]) for axis in range(3))
    if any(len(axis) < 2 for axis in axes):
        raise ValueError("offset samples must span three dimensions")
    if np.prod([len(axis) for axis in axes]) == len(points):
        order = np.lexsort(points.T[::-1])
        interpolate = RegularGridInterpolator(
            axes,
            samples[order].reshape(tuple(len(axis) for axis in axes)),
            bounds_error=False,
            fill_value=None,
        )
        nearest = None
    else:
        interpolate = LinearNDInterpolator(points, samples)
        nearest = NearestNDInterpolator(points, samples)

    def field(x, y, z):
        coordinates = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
        )
        query = np.column_stack([axis.ravel() for axis in coordinates])
        result = np.asarray(interpolate(query))
        if nearest is not None:
            outside = np.isnan(result)
            result[outside] = nearest(query[outside])
        return result.reshape(coordinates[0].shape)

    return field
