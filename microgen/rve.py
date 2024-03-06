"""
Representative Volume Element (RVE) or Representative Elementary Volume (REV)
"""

from __future__ import annotations

import cadquery as cq
import numpy as np


class Rve:
    """
    :param dim_x: X dimension of the RVE
    :param dim_y: Y dimension of the RVE
    :param dim_z: Z dimension of the RVE
    :param center: center of the RVE
    """

    def __init__(
        self,
        dim_x: float = 1,
        dim_y: float = 1,
        dim_z: float = 1,
        center: np.ndarray | tuple = (0, 0, 0),
    ) -> None:
        self.center = center
        self.dim_x = dim_x
        self.dim_y = dim_y
        self.dim_z = dim_z
        self.x_min = center[0] - 0.5 * dim_x
        self.x_max = center[0] + 0.5 * dim_x
        self.y_min = center[1] - 0.5 * dim_y
        self.y_max = center[1] + 0.5 * dim_y
        self.z_min = center[2] - 0.5 * dim_z
        self.z_max = center[2] + 0.5 * dim_z
        self.dx = abs(self.x_max - self.x_min)
        self.dy = abs(self.y_max - self.y_min)
        self.dz = abs(self.z_max - self.z_min)
        self.box = (
            cq.Workplane()
            .box(self.dim_x, self.dim_y, self.dim_z)
            .translate((self.center[0], self.center[1], self.center[2]))
        )
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
        return Rve(dim_x, dim_y, dim_z, center)
