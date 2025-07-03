"""Basic Geometry.

====================================================
Basic Geometry (:mod:`microgen.shape.shape`)
====================================================
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from scipy.spatial.transform import Rotation

if TYPE_CHECKING:
    import cadquery as cq
    import pyvista as pv

    from microgen.shape import KwargsGenerateType, Vector3DType


class ImplicitOperationUnavailableForExplicitShapeError(Exception):
    """Raised when an implicit operation is called on an explicit shape."""


class Shape(ABC):
    """Shape class to manage shapes.

    :param shape: name of the shape
    :param center: center
    :param orientation: orientation
    """

    def __init__(
        self: Shape,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType | Rotation = (0, 0, 0),
    ) -> None:
        """Initialize the shape."""
        self.center = center
        self.orientation = (
            orientation
            if isinstance(orientation, Rotation)
            else Rotation.from_euler("ZXZ", orientation, degrees=True)
        )

    @abstractmethod
    def fillet_shape(self: Shape, radius: float) -> Shape:
        """Apply a fillet to the shape.

        :param radius: radius of the fillet
        :return: Shape with fillet applied
        """
        raise NotImplementedError

    @abstractmethod
    def round_shape(self: Shape, radius: float) -> Shape:
        """Round the convex edges of the shape.

        :param radius: radius of the rounding
        :return: Shape with rounded edges
        """
        raise NotImplementedError

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


class ExplicitShape(Shape):
    """Explicit shape class to manage shapes with explicit methods.

    This class is used to define shapes that have explicit methods for
    generating the CAD shape and VTK mesh.
    """

    def fillet_shape(self: Shape, radius: float) -> ExplicitShape:
        raise ImplicitOperationUnavailableForExplicitShapeError(
            "Fillet operation is not available for explicit shapes. "
        )

    def round_shape(self: Shape, radius: float) -> ExplicitShape:
        raise ImplicitOperationUnavailableForExplicitShapeError(
            "Round operation is not available for explicit shapes. "
        )


class ImplicitShape(Shape):
    """Implicit shape class to manage shapes with implicit methods.

    This class is used to define shapes that have implicit methods for
    generating the CAD shape and VTK mesh.
    """

    def fillet_shape(self: Shape, radius: float) -> ImplicitShape: ...

    def round_shape(self: Shape, radius: float) -> ImplicitShape: ...
