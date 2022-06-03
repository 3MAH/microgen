import cadquery as cq
import pyvista as pv
import numpy as np

from ..operations import rotateEuler
from ..pvoperations import rotatePvEuler

# ----------ExtrudedPolygon-----------------------------------------------------------------------------------------#


class ExtrudedPolygon:
    """
    Class to generate an extruded polygon with a given list of points and a thickness
    """
    def __init__(
        self,
        center: np.ndarray,
        angle: np.ndarray,
        listCorners: list,
        height: float,
        number: int,
    ) -> None:
        self.center = center
        self.angle = angle
        self.listCorners = listCorners
        self.height = height
        self.number = number
        self.name_part = "extrudedpolygon" + str(self.number)

    def createExtrudedpolygon(self) -> cq.Workplane:
        poly = (
            cq.Workplane("YZ")
            .polyline(self.listCorners)
            .close()
            .extrude(self.height)
            .translate(
                (self.center[0] - self.height / 2.0, self.center[1], self.center[2])
            )
        )
        poly = rotateEuler(
            poly, self.center, self.angle[0], self.angle[1], self.angle[2]
        )
        return poly

    def createPvExtrudedpolygon(self, capping=True) -> pv.PolyData:
        vertices = []
        for corner in self.listCorners:
            vertices.append([self.center[0] - 0.5 * self.height,
                             self.center[1] + corner[0],
                             self.center[2] + corner[1]])
        faces = np.arange(len(vertices))
        faces = np.insert(faces, 0, len(vertices))

        poly = pv.PolyData(vertices, faces)
        poly = poly.extrude([self.height, 0, 0], capping=capping)

        poly = rotatePvEuler(
            poly, self.center, self.angle[0], self.angle[1], self.angle[2]
        )
        return poly
