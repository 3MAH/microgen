"""Ellipsoid.

=============================================
Ellipsoid (:mod:`microgen.shape.ellipsoid`)
=============================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pyvista as pv

from microgen.operations import rotate

from .shape import Shape

if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, Vector3DType


class Ellipsoid(Shape):
    """Class to generate an ellipsoid.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Ellipsoid().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Ellipsoid,
        radii: tuple[float, float, float] = (1, 0.5, 0.25),
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the ellipsoid."""
        super().__init__(**kwargs)
        self.radii = radii

    def generate(self: Ellipsoid, **_: KwargsGenerateType) -> CadShape:
        """Generate an ellipsoid CAD shape (OCCT).  Requires the ``[cad]`` extra."""
        from microgen.cad import make_ellipsoid

        shape = make_ellipsoid(radii=self.radii, center=self.center)
        return rotate(shape, self.center, self.orientation)

    def generate_vtk(self: Ellipsoid, **_: KwargsGenerateType) -> pv.PolyData:
        """Generate an ellipsoid VTK polydta using the given parameters."""
        transform_matrix = np.array(
            [
                [self.radii[0], 0, 0, self.center[0]],
                [0, self.radii[1], 0, self.center[1]],
                [0, 0, self.radii[2], self.center[2]],
                [0, 0, 0, 1],
            ],
        )
        sphere = pv.Sphere(radius=1)
        ellipsoid = sphere.transform(transform_matrix, inplace=False)
        return rotate(ellipsoid, self.center, self.orientation)
