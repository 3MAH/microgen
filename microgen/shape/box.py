"""Box.

===============================
Box (:mod:`microgen.shape.box`)
===============================
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import cadquery as cq
import pyvista as pv

from microgen.operations import rotateEuler, rotatePvEuler

from .shape import Shape

if TYPE_CHECKING:
    from microgen.shape import KwargsGenerateType, Vector3DType


class Box(Shape):
    """Class to generate a box.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Box().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Box,
        dim: tuple[float, float, float] = (1, 1, 1),
        dim_x: float | None = None,
        dim_y: float | None = None,
        dim_z: float | None = None,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the box."""
        super().__init__(**kwargs)
        if dim_x is not None or dim_y is not None or dim_z is not None:
            warnings.warn(
                "The 'dim_x', 'dim_y', and 'dim_z' parameters are deprecated. \
                    Use 'dim' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if dim_x is None:
                dim_x = dim[0]
            if dim_y is None:
                dim_y = dim[1]
            if dim_z is None:
                dim_z = dim[2]
            self.dim = (dim_x, dim_y, dim_z)
        self.dim = dim

    def generate(self: Box, **_: KwargsGenerateType) -> cq.Shape:
        """Generate a box CAD shape using the given parameters."""
        box = cq.Workplane().box(*self.dim).translate(self.center)
        box = rotateEuler(
            box,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )
        return cq.Shape(box.val().wrapped)

    def generate_vtk(
        self: Box,
        level: int = 0,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate a box VTK shape using the given parameters."""
        box = pv.Box(
            bounds=(
                self.center[0] - 0.5 * self.dim[0],
                self.center[0] + 0.5 * self.dim[0],
                self.center[1] - 0.5 * self.dim[1],
                self.center[1] + 0.5 * self.dim[1],
                self.center[2] - 0.5 * self.dim[2],
                self.center[2] + 0.5 * self.dim[2],
            ),
            level=level,
            quads=True,
        )
        return rotatePvEuler(
            box,
            self.center,
            self.orientation[0],
            self.orientation[1],
            self.orientation[2],
        )

    def generateVtk(self: Box, **kwargs: KwargsGenerateType) -> pv.PolyData:  # noqa: N802
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(**kwargs)
