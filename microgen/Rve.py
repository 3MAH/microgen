# -*- coding: utf8 -*-
import cadquery as cq


class Rve:
    def __init__(self, dim_x, dim_y, dim_z, size_mesh=0.03):
        """ DESCRIPTION
        
        Representative Volume Element (RVE) or Representative Elementary Volume (REV)

        Parameters
        ----------
        dim_x : float
            X dimension of the RVE
        dim_y : float
            Y dimension of the RVE
        dim_z : float
            Z dimension of the RVE
        size_mesh : float
            useless
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
