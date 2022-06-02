"""
===============================
Box (:mod:`microgen.shape.box`)
===============================
"""
import cadquery as cq
import numpy as np

from ..operations import rotateEuler


class Box:
    """
    Class to generate a box
    """
    def __init__(
        self,
        center: np.ndarray,
        orientation: np.ndarray,
        dim_x: float,
        dim_y: float,
        dim_z: float,
        number: int = 0,
    ) -> None:
        self.center = center
        self.orientation = orientation
        self.dim_x = dim_x
        self.dim_y = dim_y
        self.dim_z = dim_z
        self.number = number
        self.name_part = "box" + str(self.number)

    def createBox(self) -> cq.Workplane:
        box = (
            cq.Workplane()
            .box(self.dim_x, self.dim_y, self.dim_z)
            .translate((self.center[0], self.center[1], self.center[2]))
        )
        box = rotateEuler(box, self.center, self.orientation[0], self.orientation[1], self.orientation[2])
        return box
