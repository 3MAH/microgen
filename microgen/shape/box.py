"""Box.

===============================
Box (:mod:`microgen.shape.box`)
===============================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pyvista as pv

from microgen.operations import rotate

from .shape import Shape

if TYPE_CHECKING:
    from microgen.cad import CadShape
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
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the box."""
        super().__init__(**kwargs)
        self.dim = dim

    def generate(self: Box, **_: KwargsGenerateType) -> CadShape:
        """Generate a box CAD shape (OCCT).  Requires the ``[cad]`` extra."""
        from microgen.cad import make_box

        shape = make_box(self.dim, self.center)
        return rotate(shape, self.center, self.orientation)

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
        return rotate(box, self.center, self.orientation)
