"""Basic Geometry.

====================================================
Basic Geometry (:mod:`microgen.shape.shape`)
====================================================
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import cadquery as cq
    import pyvista as pv

    from microgen.shape import KwargsGenerateType, Vector3DType


class Shape(ABC):
    """Shape class to manage shapes.

    :param shape: name of the shape
    :param center: center
    :param orientation: orientation
    """

    def __init__(
        self: Shape,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType = (0, 0, 0),
    ) -> None:
        """Initialize the shape."""
        self.center = center
        self.orientation = orientation

    @abstractmethod
    def generate(self: Shape, **_: KwargsGenerateType) -> cq.Shape:
        """Generate the CAD shape.

        :return: cq.Shape
        """
        raise NotImplementedError

    @abstractmethod
    def generate_vtk(self: Shape, **_: KwargsGenerateType) -> pv.PolyData:
        """Generate the vtk mesh of the shape.

        :return: pv.PolyData
        """
        raise NotImplementedError

    @abstractmethod
    def generateVtk(self: Shape, **_: KwargsGenerateType) -> pv.PolyData:  # noqa: N802
        """Deprecated. Use generate_vtk instead."""  # noqa: D401
        return self.generate_vtk()
