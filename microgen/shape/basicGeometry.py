"""
====================================================
Basic Geometry (:mod:`microgen.shape.basicGeometry`)
====================================================
"""

from __future__ import annotations

from typing import Tuple

import cadquery as cq


class BasicGeometry:
    """
    BasicGeometry class to manage shapes

    :param shape: name of the shape
    :param center: center
    :param orientation: orientation
    """

    numInstances = 0

    def __init__(
        self,
        shape: str,
        center: Tuple[float, float, float] = (0, 0, 0),
        orientation: Tuple[float, float, float] = (0, 0, 0),
    ) -> None:
        self.number = self.numInstances
        self.shape = shape
        self.center = center
        self.orientation = orientation
        self.name = self.shape + "_" + str(self.number)

        self.geometry: cq.Shape | None = None
        BasicGeometry.numInstances += 1
