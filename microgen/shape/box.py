import cadquery as cq
import numpy as np

from ..operations import rotateEuler

# ----------BOX-----------------------------------------------------------------------------------------#
# MB 03/12/2021


class Box:
    def __init__(
        self,
        center: np.ndarray,
        angle: np.ndarray,
        dim_x: float,
        dim_y: float,
        dim_z: float,
        number: int,
    ) -> None:
        self.center = center
        self.angle = angle
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
        box = rotateEuler(box, self.center, self.angle[0], self.angle[1], self.angle[2])
        return box
