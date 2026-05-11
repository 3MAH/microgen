"""Capsule.

=======================================
Capsule (:mod:`microgen.shape.capsule`)
=======================================
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


class Capsule(Shape):
    """Class to generate a capsule (cylinder with hemispherical ends).

    Canonical frame: capsule axis along ``+x``, segment from
    ``(-height/2, 0, 0)`` to ``(+height/2, 0, 0)``, radius ``radius``.
    The implicit field is the point-to-segment distance minus radius
    (Inigo Quilez capsule SDF), evaluated in the local frame so
    ``center`` and ``orientation`` transform the field correctly. Set on
    every instance so capsules compose via ``|`` / ``&`` / ``-`` and
    stay usable without the ``[cad]`` extra (only :meth:`generate`
    requires CAD).

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Capsule().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(
        self: Capsule,
        height: float = 1,
        radius: float = 0.5,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the capsule."""
        super().__init__(**kwargs)
        self.height = height
        self.radius = radius
        self._setup_frep_field()

    def _setup_frep_field(self: Capsule) -> None:
        """Bake the capsule SDF and AABB onto ``_func`` / ``_bounds``."""
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
            t = np.clip(lx, -half_h, half_h)
            dx = lx - t
            return np.sqrt(dx**2 + ly**2 + lz**2) - r

        # World-space AABB: rotate the 8 corners of the canonical [(-h/2-r, h/2+r) x (-r, r) x (-r, r)] box.
        rot = self.orientation.as_matrix()
        corners = np.array(
            list(itertools.product([-half_h - r, half_h + r], [-r, r], [-r, r])),
        )
        rotated = corners @ rot.T
        margin = (half_h + r) * 0.1
        self._func = _field
        self._bounds = (
            cx + float(rotated[:, 0].min()) - margin,
            cx + float(rotated[:, 0].max()) + margin,
            cy + float(rotated[:, 1].min()) - margin,
            cy + float(rotated[:, 1].max()) + margin,
            cz + float(rotated[:, 2].min()) - margin,
            cz + float(rotated[:, 2].max()) + margin,
        )

    def generate(self: Capsule, **_: KwargsGenerateType) -> CadShape:
        """Generate a capsule CAD shape (OCCT).  Requires the ``[cad]`` extra."""
        from microgen.cad import make_capsule

        shape = make_capsule(
            radius=self.radius,
            height=self.height,
            center=self.center,
        )
        return rotate(shape, self.center, self.orientation)

    def generate_vtk(
        self: Capsule,
        resolution: int = 100,
        theta_resolution: int = 50,
        phi_resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Generate a capsule VTK shape using the given parameters."""
        cylinder = pv.Cylinder(
            center=self.center,
            radius=self.radius,
            height=self.height,
            resolution=resolution,
            capping=True,
        ).triangulate()
        sphere_left = pv.Sphere(
            radius=self.radius,
            center=(self.center[0] - self.height / 2, self.center[1], self.center[2]),
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        ).triangulate()
        sphere_right = pv.Sphere(
            radius=self.radius,
            center=(self.center[0] + self.height / 2, self.center[1], self.center[2]),
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        ).triangulate()
        capsule = cylinder.boolean_union(sphere_left).boolean_union(sphere_right)
        return rotate(capsule, self.center, self.orientation)

    def generateVtk(  # noqa: N802
        self: Capsule,
        resolution: int = 100,
        theta_resolution: int = 50,
        phi_resolution: int = 50,
        **_: KwargsGenerateType,
    ) -> pv.PolyData:
        """Deprecated. Use :meth:`generate_vtk` instead."""  # noqa: D401
        return self.generate_vtk(
            resolution=resolution,
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        )
