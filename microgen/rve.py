"""
Representative Volume Element (RVE) or Representative Elementary Volume (REV)
"""

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
        dim_x: Optional[float] = None,
        dim_y: Optional[float] = None,
        dim_z: Optional[float] = None,
        center: Union[np.ndarray, Tuple, List] = (0, 0, 0),
        dim: Optional[Union[float, np.ndarray, Tuple, List]] = 1,
    ) -> None:
        if isinstance(center, (Tuple, List)) and len(center) == _DIM:
            self.center = np.array(center)
        elif isinstance(center, np.ndarray) and center.shape == (_DIM,):
            self.center = center
        else:
            raise ValueError(f"center must be an array or Sequence of length {_DIM}")

        if (
            dim is not None
            and dim_x is not None
            and dim_y is not None
            and dim_z is not None
        ):
            raise ValueError("Either dim or dim_x, dim_y, dim_z must be specified")

        if isinstance(dim, (int, float)):
            self.dim = np.array([dim for _ in range(_DIM)])
        elif isinstance(dim, (Tuple, List)) and len(dim) == _DIM:
            self.dim = np.array(dim)
        elif isinstance(dim, np.ndarray) and dim.shape == (_DIM,):
            self.dim = dim
        elif dim is None:
            if dim_x is None or dim_y is None or dim_z is None:
                raise ValueError("dim must be specified")
            warnings.warn(
                "dim_x, dim_y, dim_z will be deprecated, use 'dim' instead.",
                DeprecationWarning,
            )
            self.dim = np.array([dim_x, dim_y, dim_z])
        else:
            raise ValueError(f"dim must be an array or Sequence of length {_DIM}")

        self.min_point = self.center - 0.5 * self.dim
        self.max_point = self.center + 0.5 * self.dim

        self.box = cq.Workplane().box(*self.dim).translate(vec=cq.Vector(*self.center))
        self.is_matrix = False
        self.matrix_number = 0

        # Deprecated attributes
        if dim_x is not None and dim_y is not None and dim_z is not None:
            self.dim_x = dim_x
            self.dim_y = dim_y
            self.dim_z = dim_z
        self.x_min = center[0] - 0.5 * self.dim[0]
        self.x_max = center[0] + 0.5 * self.dim[0]
        self.y_min = center[1] - 0.5 * self.dim[1]
        self.y_max = center[1] + 0.5 * self.dim[1]
        self.z_min = center[2] - 0.5 * self.dim[2]
        self.z_max = center[2] + 0.5 * self.dim[2]
        self.dx = abs(self.x_max - self.x_min)
        self.dy = abs(self.y_max - self.y_min)
        self.dz = abs(self.z_max - self.z_min)
