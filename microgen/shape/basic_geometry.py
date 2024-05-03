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

    from microgen.shape import KwargsGenerateType, Vector3DType


class BasicGeometry:
    """BasicGeometry class to manage shapes.

    :param shape: name of the shape
    :param center: center
    :param orientation: orientation
    """

    num_instances = 0

    def __init__(
        self: BasicGeometry,
        shape: str,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType = (0, 0, 0),
    ) -> None:
        """Initialize the shape."""
        self.number = self.num_instances
        self.shape = shape
        self.center = center
        self.orientation = orientation
        self.name = f"{self.shape}_{self.number}"

        self.geometry: cq.Shape | None = None
        BasicGeometry.num_instances += 1

    def generate(self: BasicGeometry, **_: KwargsGenerateType) -> cq.Shape:
        """Generate the CAD shape.

        :return: cq.Shape
        """
        raise NotImplementedError

    def generate_vtk(self: BasicGeometry, **_: KwargsGenerateType) -> pv.PolyData:
        """Generate the vtk mesh of the shape.

        :return: pv.PolyData
        """
        raise NotImplementedError

    def generateVtk(self: BasicGeometry, **_: KwargsGenerateType) -> pv.PolyData:  # noqa: N802
        """Deprecated. Use generate_vtk instead."""  # noqa: D401
        return self.generate_vtk()
