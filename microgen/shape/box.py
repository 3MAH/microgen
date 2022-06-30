"""
===============================
Box (:mod:`microgen.shape.box`)
===============================
"""
import cadquery as cq
import pyvista as pv

from ..operations import rotateEuler, rotatePvEuler
from .basicGeometry import BasicGeometry


class Box(BasicGeometry):
    """
    Class to generate a box
    """

    def __init__(
        self,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        dim_x: float = 1,
        dim_y: float = 1,
        dim_z: float = 1,
    ) -> None:
        super().__init__(shape="Box", center=center, orientation=orientation)
        self.dim_x = dim_x
        self.dim_y = dim_y
        self.dim_z = dim_z

    def generate(self) -> cq.Shape:
        box = (
            cq.Workplane()
            .box(self.dim_x, self.dim_y, self.dim_z)
            .translate((self.center[0], self.center[1], self.center[2]))
        )
        box = rotateEuler(
            box,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return cq.Shape(box.val().wrapped)

    def generateVtk(self, level=0, quads=True) -> pv.PolyData:
        box = pv.Box(
            bounds=(
                self.center[0] - 0.5 * self.dim_x,
                self.center[0] + 0.5 * self.dim_x,
                self.center[1] - 0.5 * self.dim_y,
                self.center[1] + 0.5 * self.dim_y,
                self.center[2] - 0.5 * self.dim_z,
                self.center[2] + 0.5 * self.dim_z,
            ),
            level=level,
            quads=quads,
        )
        box = rotatePvEuler(
            box,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return box
