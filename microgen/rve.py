"""
Representative Volume Element (RVE) or Representative Elementary Volume (REV)
"""

from __future__ import annotations

import warnings
from typing import List, Optional, Tuple, Union

import cadquery as cq
import numpy as np

_DIM = 3


class Rve:
    """
    :param dim_x: X dimension of the RVE
    :param dim_y: Y dimension of the RVE
    :param dim_z: Z dimension of the RVE
    :param center: center of the RVE
    :param dim: dimensions of the RVE
    """

    def __init__(
        self,
        center: Union[np.ndarray, Tuple, List] = (0, 0, 0),
        dim: Optional[Union[float, np.ndarray, Tuple, List]] = 1,
    ) -> None:
        if isinstance(center, (Tuple, List)) and len(center) == _DIM:
            self.center = np.array(center)
        elif isinstance(center, np.ndarray) and center.shape == (_DIM,):
            self.center = center
        else:
            raise ValueError(f"center must be an array or Sequence of length {_DIM}")

        if isinstance(dim, (int, float)):
            self.dim = np.array([dim for _ in range(_DIM)])
        elif isinstance(dim, (Tuple, List)) and len(dim) == _DIM:
            self.dim = np.array(dim)
        elif isinstance(dim, np.ndarray) and dim.shape == (_DIM,):
            self.dim = dim
        else:
            raise ValueError(f"dim must be an array or Sequence of length {_DIM}")

        self.min_point = self.center - 0.5 * self.dim
        self.max_point = self.center + 0.5 * self.dim

        self.box = cq.Workplane().box(*self.dim).translate(vec=cq.Vector(*self.center))
        self.is_matrix = False
        self.matrix_number = 0

    @classmethod
    def from_min_max(
        cls,
        x_min: float = -0.5,
        x_max: float = 0.5,
        y_min: float = -0.5,
        y_max: float = 0.5,
        z_min: float = -0.5,
        z_max: float = 0.5,
    ) -> Rve:
        """
        :param x_min: min X dimension of the RVE
        :param x_max: max X dimension of the RVE
        :param x_min: min Y dimension of the RVE
        :param x_max: max Y dimension of the RVE
        :param x_min: min Z dimension of the RVE
        :param x_max: max Z dimension of the RVE
        Generate a Rve from the min - max values
        """

        center = (0.5 * (x_min + x_max), 0.5 * (y_min + y_max), 0.5 * (z_min + z_max))
        dim_x = abs(x_max - x_min)
        dim_y = abs(y_max - y_min)
        dim_z = abs(z_max - z_min)
        return Rve(dim=(dim_x, dim_y, dim_z), center=center)
