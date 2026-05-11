"""Representative Volume Element (RVE).

The ``Rve.box`` attribute is a :class:`~microgen.cad.CadShape` wrapping an
OCCT box; it is built lazily on first access and requires the ``[cad]``
extra (``cadquery-ocp-novtk``).
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np

_DIM = 3

if TYPE_CHECKING:
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

        self.is_matrix = False
        self.matrix_number = 0

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
        x_min: float = -0.5,
        x_max: float = 0.5,
        y_min: float = -0.5,
        y_max: float = 0.5,
        z_min: float = -0.5,
        z_max: float = 0.5,
    ) -> Rve:
        """Generate a Rve from the min - max values.

        :param x_min: min X dimension of the RVE
        :param x_max: max X dimension of the RVE
        :param x_min: min Y dimension of the RVE
        :param x_max: max Y dimension of the RVE
        :param x_min: min Z dimension of the RVE
        :param x_max: max Z dimension of the RVE
        """
        center = (0.5 * (x_min + x_max), 0.5 * (y_min + y_max), 0.5 * (z_min + z_max))
        dim_x = abs(x_max - x_min)
        dim_y = abs(y_max - y_min)
        dim_z = abs(z_max - z_min)
        return Rve(dim=(dim_x, dim_y, dim_z), center=center)
