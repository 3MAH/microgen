"""Ellipsoid.

=============================================
Ellipsoid (:mod:`microgen.shape.ellipsoid`)
=============================================
"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pyvista as pv

from microgen.operations import rotate

from .shape import Shape

if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, Vector3DType


class Ellipsoid(Shape):
    """Class to generate an ellipsoid.

    The implicit field is the canonical ellipsoid scalar
    ``sqrt((x/rx)^2 + (y/ry)^2 + (z/rz)^2) - 1`` (approximate SDF with
    the correct zero-level set), evaluated in the local frame so
    ``center`` and ``orientation`` transform the field correctly. Set
    on every instance so ellipsoids compose via ``|`` / ``&`` / ``-``
    and stay usable without the ``[cad]`` extra (only :meth:`generate_cad`
    requires CAD).

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Ellipsoid().generate_surface_mesh()
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
        self._setup_frep_field()

    def _setup_frep_field(self: Ellipsoid) -> None:
        """Bake the ellipsoid field and AABB onto ``_func`` / ``_bounds``."""
        cx, cy, cz = (float(c) for c in self.center)
        rx, ry, rz = (float(r) for r in self.radii)
        rot_inv = self.orientation.inv().as_matrix()

        def _field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            px = x - cx
            py = y - cy
            pz = z - cz
            lx = rot_inv[0, 0] * px + rot_inv[0, 1] * py + rot_inv[0, 2] * pz
            ly = rot_inv[1, 0] * px + rot_inv[1, 1] * py + rot_inv[1, 2] * pz
            lz = rot_inv[2, 0] * px + rot_inv[2, 1] * py + rot_inv[2, 2] * pz
            return np.sqrt((lx / rx) ** 2 + (ly / ry) ** 2 + (lz / rz) ** 2) - 1.0

        rot = self.orientation.as_matrix()
        corners = np.array(list(itertools.product([-rx, rx], [-ry, ry], [-rz, rz])))
        rotated = corners @ rot.T
        margin = max(rx, ry, rz) * 0.1
        self._func = _field
        self._bounds = (
            cx + float(rotated[:, 0].min()) - margin,
            cx + float(rotated[:, 0].max()) + margin,
            cy + float(rotated[:, 1].min()) - margin,
            cy + float(rotated[:, 1].max()) + margin,
            cz + float(rotated[:, 2].min()) - margin,
            cz + float(rotated[:, 2].max()) + margin,
        )

    def generate_cad(self: Ellipsoid, **_: KwargsGenerateType) -> CadShape:
        """Generate an ellipsoid CAD shape (OCCT).  Requires the ``[cad]`` extra."""
        from microgen.cad import make_ellipsoid

        shape = make_ellipsoid(radii=self.radii, center=self.center)
        return rotate(shape, self.center, self.orientation)

    def generate_surface_mesh(self: Ellipsoid, **_: KwargsGenerateType) -> pv.PolyData:
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
