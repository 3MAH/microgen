"""Representative Volume Element (RVE).

The ``Rve.box`` attribute is a :class:`~microgen.cad.CadShape` wrapping an
OCCT box; it is built lazily on first access and requires the ``[cad]``
extra (``cadquery-ocp-novtk``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

_DIM = 3

if TYPE_CHECKING:
    from collections.abc import Sequence

    from .cad import CadShape

    Vector3DType = tuple[float, float, float] | Sequence[float]


class Rve:
    """Representative Volume Element (RVE).

    :param center: center of the RVE
    :param dim: dimensions of the RVE
    """

    def __init__(
        self: Rve,
        center: Vector3DType = (0, 0, 0),
        dim: float | Vector3DType = 1,
    ) -> None:
        """Initialize the RVE."""
        if isinstance(center, (tuple, list)) and len(center) == _DIM:
            self.center = np.array(center)
        elif isinstance(center, np.ndarray) and center.shape == (_DIM,):
            self.center = center
        else:
            err_msg = f"center must be an array or Sequence of length {_DIM}"
            raise ValueError(err_msg)

        if isinstance(dim, (int, float)):
            self.dim = np.array([dim for _ in range(_DIM)])
        elif isinstance(dim, (tuple, list)) and len(dim) == _DIM:
            self.dim = np.array(dim)
        elif isinstance(dim, np.ndarray) and dim.shape == (_DIM,):
            self.dim = dim
        else:
            err_msg = f"dim must be an array or Sequence of length {_DIM}"
            raise ValueError(err_msg)

        if np.any(self.dim <= 0):
            err_msg = f"dimensions of the RVE must be greater than 0, got {self.dim}"
            raise ValueError(err_msg)

        self.min_point = self.center - 0.5 * self.dim
        self.max_point = self.center + 0.5 * self.dim

        self._cached_box: CadShape | None = None

    @property
    def box(self) -> CadShape:
        """Return a :class:`~microgen.cad.CadShape` box of the RVE (cached).

        Requires the ``[cad]`` extra.
        """
        if self._cached_box is None:
            from .cad import make_box  # noqa: PLC0415

            self._cached_box = make_box(tuple(self.dim), tuple(self.center))
        return self._cached_box

    @box.setter
    def box(self, value: CadShape) -> None:
        self._cached_box = value

    @classmethod
    def from_min_max(
        cls: type[Rve],
        min_point: Vector3DType = (-0.5, -0.5, -0.5),
        max_point: Vector3DType = (0.5, 0.5, 0.5),
    ) -> Rve:
        """Generate a Rve from min and max corner points.

        :param min_point: ``(x_min, y_min, z_min)`` corner of the RVE
        :param max_point: ``(x_max, y_max, z_max)`` corner of the RVE
        """
        lo = np.asarray(min_point, dtype=float)
        hi = np.asarray(max_point, dtype=float)
        return cls(center=tuple(0.5 * (lo + hi)), dim=tuple(np.abs(hi - lo)))
