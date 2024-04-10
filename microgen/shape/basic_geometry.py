"""Basic Geometry.

====================================================
Basic Geometry (:mod:`microgen.shape.basic_geometry`)
====================================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import cadquery as cq
    import pyvista as pv


class BasicGeometry:
    """BasicGeometry class to manage shapes.

    :param shape: name of the shape
    :param center: center
    :param orientation: orientation
    """

    num_instances = 0

    def __init__(
        self,
        shape: str,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
    ) -> None:
        """Initialize the shape."""
        self.number = self.num_instances
        self.shape = shape
        self.center = center
        self.orientation = orientation
        self.name = f"{self.shape}_{self.number}"

        self.geometry: cq.Shape | None = None
        BasicGeometry.num_instances += 1

    def generate(self, **_) -> cq.Shape:
        """Generate the CAD shape.

        :return: cq.Shape
        """
        raise NotImplementedError

    def generate_vtk(self, **_) -> pv.PolyData:
        """Generate the vtk mesh of the shape.

        :return: pv.PolyData
        """
        raise NotImplementedError

    def generateVtk(self, **_) -> pv.PolyData:  # noqa: N802
        """Deprecated. Use generate_vtk instead."""  # noqa: D401
        return self.generate_vtk(**_)
