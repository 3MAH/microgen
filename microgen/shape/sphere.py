"""Sphere.

=====================================
Sphere (:mod:`microgen.shape.sphere`)
=====================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import pyvista as pv

from .shape import Shape

if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, Vector3DType


class Sphere(Shape):
    """Class to generate a sphere.

    The implicit field ``f(p) = ||p - center|| - radius`` (a true SDF,
    negative inside) is set on every instance so spheres compose with
    other shapes through ``|`` / ``&`` / ``-`` and stay usable when the
    ``[cad]`` extra is not installed (only :meth:`generate` requires CAD).

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Sphere().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Sphere,
        radius: float = 1,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the sphere."""
        super().__init__(**kwargs)
        self.radius = radius
        self._setup_frep_field()

    def _setup_frep_field(self: Sphere) -> None:
        """Bake the sphere SDF and AABB onto ``_func`` / ``_bounds``."""
        cx, cy, cz = (float(c) for c in self.center)
        r = float(self.radius)
        margin = r * 1.1

        def _field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            return np.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) - r

        self._func = _field
        self._bounds = (
            cx - margin,
            cx + margin,
            cy - margin,
            cy + margin,
            cz - margin,
            cz + margin,
        )

    def generate(self: Sphere, **_: KwargsGenerateType) -> CadShape:
        """Generate a sphere CAD shape (OCCT).  Requires the ``[cad]`` extra."""
        from microgen.cad import make_sphere

        return make_sphere(self.radius, self.center)

    def generate_vtk(
        self: Sphere,
        theta_resolution: int = 50,
        phi_resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate a sphere VTK shape using the given parameters."""
        return pv.Sphere(
            radius=self.radius,
            center=tuple(self.center),
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        )
