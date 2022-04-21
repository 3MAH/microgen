# -*- coding: utf8 -*-
import cadquery as cq


class Rve:
    def __init__(self, dim_x, dim_y, dim_z, size_mesh):
    """

    Parameters
    ----------
    dim_x : TYPE
        DESCRIPTION
    dim_y : TYPE
        DESCRIPTION
    dim_z : TYPE
        DESCRIPTION
    size_mesh : TYPE
        DESCRIPTION
    """
        self.dim_x = dim_x
        self.dim_y = dim_y
        self.dim_z = dim_z
        self.x_min = 0.0
        self.x_max = dim_x
        self.y_min = 0.0
        self.y_max = dim_y
        self.z_min = 0.0
        self.z_max = dim_z
        self.size_mesh = size_mesh
        self.dx = abs(self.x_max - self.x_min)
        self.dy = abs(self.y_max - self.y_min)
        self.dz = abs(self.z_max - self.z_min)
        self.Box = cq.Workplane("XY").box(
            self.dx, self.dy, self.dz, centered=(False, False, False)
        )
        self.is_matrix = False
        self.matrix_number = 0
