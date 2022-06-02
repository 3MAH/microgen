"""
========================================================
Extruded Polygon (:mod:`microgen.shape.extrudedPolygon`)
========================================================
"""
import cadquery as cq
import numpy as np

from typing import Sequence, Tuple

from ..operations import rotateEuler


class ExtrudedPolygon:
    """
    Class to generate an extruded polygon with a given list of points and a thickness
    """
    def __init__(
        self,
        center: np.ndarray,
        orientation: np.ndarray,
        listCorners: Sequence[Tuple[float, float]],
        height: float,
        number: int = 0,
    ) -> None:
        self.center = center
        self.orientation = orientation
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
            poly, self.center, self.orientation[0], self.orientation[1], self.orientation[2]
        )
        return poly
