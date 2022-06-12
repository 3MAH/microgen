"""
========================================================
Extruded Polygon (:mod:`microgen.shape.extrudedPolygon`)
========================================================
"""
from typing import Sequence, Tuple

import cadquery as cq
import numpy as np

from ..operations import rotateEuler
from .basicGeometry import BasicGeometry


class ExtrudedPolygon(BasicGeometry):
    """
    Class to generate an extruded polygon with a given list of points and a thickness
    """

    def __init__(
        self,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        listCorners: Sequence[Tuple[float, float]] = [
            (1, 0),
            (0.5, 0.5 * np.sqrt(3)),
            (-0.5, 0.5 * np.sqrt(3)),
            (-1, 0),
            (-0.5, -0.5 * np.sqrt(3)),
            (0.5, -0.5 * np.sqrt(3)),
            (1, 0),
        ],  # hexagon
        height: float = 1,
    ) -> None:
        super().__init__(
            shape="ExtrudedPolygon", center=center, orientation=orientation
        )
        self.listCorners = listCorners
        self.height = height

    def generate(self) -> cq.Shape:
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
            poly,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return cq.Shape(poly.val().wrapped)
