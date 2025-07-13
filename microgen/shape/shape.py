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


class ShellCreationError(Exception):
    """Error raised when the shell creation fails."""

    def __init__(self: ShellCreationError, message: str) -> None:
        """Initialize the ShellCreationError."""
        super().__init__(message)


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
        resolution: int = 20,
    ) -> None:
        """Initialize the implicit shape."""
        super().__init__(center, orientation)
        self.resolution = resolution
        self.grid: pv.StructuredGrid
        self._surface: pv.PolyData = None

    @property
    @abstractmethod
    def surface_function(self: ImplicitShape) -> Field:
        """Get the surface function of the implicit shape.

        :return: Field representing the surface function
        """
        raise NotImplementedError

    @property
    def surface(self: ImplicitShape) -> pv.PolyData:
        """Returns isosurface f(x, y, z) = 0."""
        if self._surface is not None:
            return self._surface

        self._surface = self.grid.contour(
            isosurfaces=[0.0],
            scalars="surface",
        ).triangulate()
        return self._surface

    def _create_grid(
        self: ImplicitShape,
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> pv.StructuredGrid:
        """Return the structured cartesian grid of the ImplicitShape."""
        grid = pv.StructuredGrid(x, y, z)
        grid["coords"] = np.c_[
            x.ravel(order="F"),
            y.ravel(order="F"),
            z.ravel(order="F"),
        ]
        return grid

    def _create_shell(self: ImplicitShape, mesh: pv.PolyData) -> cq.Shell:
        if not mesh.is_all_triangles:
            mesh.triangulate(inplace=True)  # useless ?
        triangles = mesh.faces.reshape(-1, 4)[:, 1:]
        triangles = np.c_[triangles, triangles[:, 0]]

        faces = []
        for tri in triangles:
            lines = [
                cq.Edge.makeLine(
                    cq.Vector(*mesh.points[start]),
                    cq.Vector(*mesh.points[end]),
                )
                for start, end in zip(tri[:], tri[1:])
            ]

            wire = cq.Wire.assembleEdges(lines)
            faces.append(cq.Face.makeFromWires(wire))

        try:
            shell = cq.Shell.makeShell(faces)
        except ValueError as err:
            err_msg = "Failed to create the shell, \
                try to increase the resolution or the smoothing."
            raise ShellCreationError(err_msg) from err
        return shell

    def _create_surface(
        self: ImplicitShape,
        isovalue: float | npt.NDArray[np.float64] = 0.0,
        smoothing: int = 0,
    ) -> cq.Shell:
        """Create an implicit surface for the given isovalue."""
        if isinstance(isovalue, (int, float)):
            scalars = self.grid["surface"] - isovalue
        elif isinstance(isovalue, np.ndarray):
            scalars = self.grid["surface"] - isovalue.ravel(order="F")

        mesh = self.grid.contour(isosurfaces=[0.0], scalars=scalars)
        mesh.smooth(n_iter=smoothing, feature_smoothing=True, inplace=True)
        mesh.clean(inplace=True)

        return self._create_shell(mesh=mesh)

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
