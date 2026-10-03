"""Bound arbitrary implicit solids explicitly and place them in a build frame."""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import TYPE_CHECKING, Literal

import numpy as np

from .types import FloatArray

if TYPE_CHECKING:
    from microgen.shape._types import BoundsType, Field
    from microgen.shape.shape import Shape


@dataclass(frozen=True, eq=False)
class BoundedSolid:
    """Negative-inside callable intersected with local bounds and transformed.

    All coordinates use millimeters. ``build_transform`` maps local points to
    build points. Shape renderer placement is deliberately not inferred.
    ``kind`` describes the supplied field, not the clipped/transformed result.
    """

    func: Field
    bounds: BoundsType
    build_transform: FloatArray = field(default_factory=lambda: np.eye(4))
    kind: Literal["levelset", "distance_estimate", "sdf"] = "levelset"
    _inverse: FloatArray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        bounds = np.asarray(self.bounds, dtype=float)
        if bounds.shape != (6,) or not np.isfinite(bounds).all():
            raise ValueError("Bounds must contain six finite coordinates")
        if np.any(bounds[1::2] <= bounds[::2]):
            raise ValueError("Each upper bound must exceed its lower bound")
        if not callable(self.func) or self.kind not in {
            "levelset",
            "distance_estimate",
            "sdf",
        }:
            raise ValueError("Require a callable field and a supported field kind")
        transform = np.array(self.build_transform, dtype=float, copy=True)
        if transform.shape != (4, 4) or not np.isfinite(transform).all():
            raise ValueError("Build transform must be a finite 4x4 affine matrix")
        if not np.array_equal(transform[3], [0, 0, 0, 1]):
            raise ValueError("Build transform must be affine")
        try:
            inverse = np.linalg.inv(transform)
        except np.linalg.LinAlgError as error:
            raise ValueError("Build transform must be invertible") from error
        transform.setflags(write=False)
        inverse.setflags(write=False)
        object.__setattr__(self, "bounds", tuple(float(value) for value in bounds))
        object.__setattr__(self, "build_transform", transform)
        object.__setattr__(self, "_inverse", inverse)

    @classmethod
    def from_shape(
        cls, shape: Shape, *, build_transform: FloatArray | None = None
    ) -> BoundedSolid:
        """Adapt a Shape's local field, requiring explicit finite clipping bounds."""
        if shape.bounds is None:
            raise ValueError("Slicing a Shape requires explicit bounds")
        return cls(
            shape.require_func(),
            shape.bounds,
            np.eye(4) if build_transform is None else build_transform,
        )

    @property
    def build_bounds(self) -> BoundsType:
        """Axis-aligned bounds of the transformed clipping box."""
        corners = np.array(list(product(*[self.bounds[i : i + 2] for i in (0, 2, 4)])))
        transformed = (
            corners @ self.build_transform[:3, :3].T + self.build_transform[:3, 3]
        )
        low, high = transformed.min(axis=0), transformed.max(axis=0)
        return (
            float(low[0]),
            float(high[0]),
            float(low[1]),
            float(high[1]),
            float(low[2]),
            float(high[2]),
        )

    def evaluate(self, x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:
        """Evaluate clipped geometry at broadcastable build coordinates.

        Field values outside the local box are never evaluated. No physical
        distance guarantee is made for these sign-preserving clipping values.
        """
        x, y, z = np.broadcast_arrays(x, y, z)
        points = np.stack((x.ravel(), y.ravel(), z.ravel()), axis=1)
        if not np.isfinite(points).all():
            raise ValueError("Sampling coordinates must be finite")
        local = points @ self._inverse[:3, :3].T + self._inverse[:3, 3]
        bounds = np.asarray(self.bounds)
        center = (bounds[::2] + bounds[1::2]) / 2
        half = (bounds[1::2] - bounds[::2]) / 2
        envelope = np.max(np.abs(local - center) - half, axis=1)
        values = envelope.copy()
        within = envelope <= 0
        if within.any():
            sample = local[within]
            raw = np.asarray(
                self.func(sample[:, 0], sample[:, 1], sample[:, 2]), dtype=float
            )
            raw = np.broadcast_to(raw, (len(sample),))
            if not np.isfinite(raw).all():
                raise ValueError(
                    "Implicit field returned nonfinite values inside its bounds"
                )
            values[within] = np.maximum(raw, envelope[within])
        return values.reshape(x.shape)
