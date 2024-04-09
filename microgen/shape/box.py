"""Box.

===============================
Box (:mod:`microgen.shape.box`)
===============================
"""

from __future__ import annotations

import cadquery as cq
import pyvista as pv

from microgen.operations import rotateEuler, rotatePvEuler

from .basic_geometry import BasicGeometry


class Box(BasicGeometry):
    """Class to generate a box.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Box().generateVtk()
       shape.plot(color='white')
    """

    def __init__(
        self,
        center: tuple[float, float, float] = (0, 0, 0),
        orientation: tuple[float, float, float] = (0, 0, 0),
        dim: tuple[float, float, float] = (1, 1, 1),
    ) -> None:
        """Initialize the box."""
        super().__init__(shape="Box", center=center, orientation=orientation)
        self.dim = dim

    def generate(self, **_) -> cq.Shape:
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
        self,
        level: int = 0,
        **_,
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

    def generateVtk(self, **kwargs) -> pv.PolyData:  # noqa: N802
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(**kwargs)
