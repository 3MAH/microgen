# -*- coding: utf8 -*-
import cadquery as cq


class Rve:
    def __init__(self, a, b, c, size_mesh):
        self.a = a
        self.b = b
        self.c = c
        self.x_min = 0.0
        self.x_max = a
        self.y_min = 0.0
        self.y_max = b
        self.z_min = 0.0
        self.z_max = c
        self.size_mesh = size_mesh
        self.dx = abs(self.x_max - self.x_min)
        self.dy = abs(self.y_max - self.y_min)
        self.dz = abs(self.z_max - self.z_min)
        self.Box = cq.Workplane("XY").box(
            self.dx, self.dy, self.dz, centered=(False, False, False)
        )
        self.is_matrix = False
        self.matrix_number = 0
