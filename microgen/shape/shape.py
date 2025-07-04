"""Basic Geometry.

====================================================
Basic Geometry (:mod:`microgen.shape.shape`)
====================================================
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Callable

import numpy as np
import numpy.typing as npt
from scipy.spatial.transform import Rotation

if TYPE_CHECKING:
    import cadquery as cq
    import pyvista as pv

    from microgen.shape import KwargsGenerateType, Vector3DType


Field = Callable[
    [npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]],
    npt.NDArray[np.float64],
]


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

    def __init__(
        self: ExplicitShape,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType | Rotation = (0, 0, 0),
    ) -> None:
        """Initialize the explicit shape."""
        super().__init__(center, orientation)


class ImplicitShape(Shape):
    """Implicit shape class to manage shapes with implicit methods.

    This class is used to define shapes that have implicit methods for
    generating the CAD shape and VTK mesh.
    """

    def __init__(
        self: ImplicitShape,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType | Rotation = (0, 0, 0),
    ) -> None:
        """Initialize the implicit shape."""
        super().__init__(center, orientation)

    @property
    @abstractmethod
    def surface_function(self: ImplicitShape) -> Field:
        """Get the surface function of the implicit shape.

        :return: Field representing the surface function
        """
        raise NotImplementedError

    ##TODO: maybe fillet and round should be in operations.py and return CustomImplicitShape
    def fillet_shape(self: ImplicitShape, radius: float) -> Field:
        """Create a fillet shape with the given radius.

        :param radius: radius of the fillet
        :return: Field representing the filleted surface
        """

        def filleted_field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            return self.surface_function(x, y, z) - radius

        return filleted_field

    def round_shape(self: ImplicitShape, radius: float) -> Field:
        """Create a rounded shape with the given radius.

        :param radius: radius of the rounding
        :return: ImplicitShape with rounding
        """

        def rounded_field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            return self.surface_function(x, y, z) + radius

        return rounded_field
