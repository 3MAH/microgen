"""
====================================================
Basic Geometry (:mod:`microgen.shape.basicGeometry`)
====================================================
"""

import cadquery as cq
import pyvista as pv

from typing import Union


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
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
    ) -> None:

        self.number = self.numInstances
        self.shape = shape
        self.center = center
        self.orientation = orientation
        self.name = self.shape + "_" + str(self.number)

        self.geometry = None  # type: Union[cq.Shape, None]
        BasicGeometry.numInstances += 1

    def generate(self) -> cq.Shape:
        """
        Generates the shape

        :return: cq.Shape
        """

    def generateVtk(self) -> pv.PolyData:
        """
        Generates the shape in vtk format

        :return: pv.PolyData
        """
