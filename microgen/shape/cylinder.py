"""
=========================================
Cylinder (:mod:`microgen.shape.cylinder`)
=========================================
"""
import cadquery as cq
import pyvista as pv

from ..operations import rotateEuler, rotatePvEuler

from .basicGeometry import BasicGeometry


class Cylinder(BasicGeometry):
    """
    Class to generate a cylinder
    """

    def __init__(
        self,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        height: float = 1,
        radius: float = 0.5,
    ) -> None:
        super().__init__(shape="Cylinder", center=center, orientation=orientation)
        self.radius = radius
        self.height = height

    def generate(self) -> cq.Shape:
        cylinder = (
            cq.Workplane("YZ")
            .circle(self.radius)
            .extrude(self.height)
            .translate(
                (self.center[0] - self.height / 2.0, self.center[1], self.center[2])
            )
        )
        cylinder = rotateEuler(
            cylinder,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return cq.Shape(cylinder.val().wrapped)

    def generateVtk(self, resolution=100, capping=True) -> pv.PolyData:
        cylinder = pv.Cylinder(
            center=tuple(self.center),
            direction=(1.0, 0.0, 0.0),
            radius=self.radius,
            height=self.height,
            resolution=resolution,
            capping=capping,
        )
        cylinder = rotatePvEuler(
            cylinder,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return cylinder
