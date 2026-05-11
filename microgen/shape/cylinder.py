"""Cylinder.

=========================================
Cylinder (:mod:`microgen.shape.cylinder`)
=========================================
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


class Cylinder(Shape):
    """Class to generate a cylinder.

    Canonical frame: cylinder axis along ``+x`` from
    ``(-height/2, 0, 0)`` to ``(+height/2, 0, 0)``, radius ``radius``.
    The implicit field is a capped-cylinder SDF (radial + axial
    distance composition), evaluated in the local frame so ``center``
    and ``orientation`` transform the field correctly. Set on every
    instance so cylinders compose via ``|`` / ``&`` / ``-`` and stay
    usable without the ``[cad]`` extra (only :meth:`generate` requires
    CAD).

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Cylinder().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Cylinder,
        height: float = 1,
        radius: float = 0.5,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the cylinder."""
        super().__init__(**kwargs)
        self.radius = radius
        self.height = height
        self._setup_frep_field()

    def _setup_frep_field(self: Cylinder) -> None:
        """Bake the cylinder SDF and AABB onto ``_func`` / ``_bounds``."""
        cx, cy, cz = (float(c) for c in self.center)
        h = float(self.height)
        r = float(self.radius)
        half_h = 0.5 * h
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
            d_radial = np.sqrt(ly**2 + lz**2) - r
            d_axial = np.abs(lx) - half_h
            outside = np.sqrt(
                np.maximum(d_radial, 0.0) ** 2 + np.maximum(d_axial, 0.0) ** 2,
            )
            inside = np.minimum(np.maximum(d_radial, d_axial), 0.0)
            return outside + inside

        rot = self.orientation.as_matrix()
        corners = np.array(
            list(itertools.product([-half_h, half_h], [-r, r], [-r, r])),
        )
        rotated = corners @ rot.T
        margin = max(r, half_h) * 0.1
        self._func = _field
        self._bounds = (
            cx + float(rotated[:, 0].min()) - margin,
            cx + float(rotated[:, 0].max()) + margin,
            cy + float(rotated[:, 1].min()) - margin,
            cy + float(rotated[:, 1].max()) + margin,
            cz + float(rotated[:, 2].min()) - margin,
            cz + float(rotated[:, 2].max()) + margin,
        )

    def generate(self: Cylinder, **_: KwargsGenerateType) -> CadShape:
        """Generate a cylinder CAD shape (OCCT).  Requires the ``[cad]`` extra."""
        from microgen.cad import make_cylinder

        shape = make_cylinder(
            radius=self.radius,
            height=self.height,
            center=self.center,
            axis=(1.0, 0.0, 0.0),
        )
        return rotate(shape, self.center, self.orientation)

    def generate_vtk(
        self: Cylinder,
        resolution: int = 100,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate a cylinder VTK shape using the given parameters."""
        cylinder = pv.Cylinder(
            center=tuple(self.center),
            direction=(1.0, 0.0, 0.0),
            radius=self.radius,
            height=self.height,
            resolution=resolution,
            capping=True,
        )
        return rotate(cylinder, self.center, self.orientation)
