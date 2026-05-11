"""Box.

===============================
Box (:mod:`microgen.shape.box`)
===============================
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


class Box(Shape):
    """Class to generate a box.

    The implicit field is the canonical AABB SDF
    (``max(|x|-hx, |y|-hy, |z|-hz)`` decomposed into outside/inside parts),
    evaluated in the box's local frame so ``center`` and ``orientation``
    transform the field correctly. Set on every instance so boxes compose
    via ``|`` / ``&`` / ``-`` and stay usable without the ``[cad]`` extra
    (only :meth:`generate` requires CAD).

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Box().generate_surface_mesh()
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
        self._setup_frep_field()

    def _setup_frep_field(self: Box) -> None:
        """Bake the box SDF and AABB onto ``_func`` / ``_bounds``."""
        cx, cy, cz = (float(c) for c in self.center)
        hx, hy, hz = (0.5 * float(d) for d in self.dim)
        rot_inv = self.orientation.inv().as_matrix()

        def _field(
            x: npt.NDArray[np.float64],
            y: npt.NDArray[np.float64],
            z: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            # Transform world coords -> box-local coords.
            px = x - cx
            py = y - cy
            pz = z - cz
            lx = rot_inv[0, 0] * px + rot_inv[0, 1] * py + rot_inv[0, 2] * pz
            ly = rot_inv[1, 0] * px + rot_inv[1, 1] * py + rot_inv[1, 2] * pz
            lz = rot_inv[2, 0] * px + rot_inv[2, 1] * py + rot_inv[2, 2] * pz
            qx = np.abs(lx) - hx
            qy = np.abs(ly) - hy
            qz = np.abs(lz) - hz
            outside = np.sqrt(
                np.maximum(qx, 0.0) ** 2
                + np.maximum(qy, 0.0) ** 2
                + np.maximum(qz, 0.0) ** 2,
            )
            inside = np.minimum(np.maximum(qx, np.maximum(qy, qz)), 0.0)
            return outside + inside

        # World-space AABB: take the 8 canonical-frame corners, rotate, recenter.
        rot = self.orientation.as_matrix()
        corners = np.array(list(itertools.product([-hx, hx], [-hy, hy], [-hz, hz])))
        rotated = corners @ rot.T
        margin = max(hx, hy, hz) * 0.1
        self._func = _field
        self._bounds = (
            cx + float(rotated[:, 0].min()) - margin,
            cx + float(rotated[:, 0].max()) + margin,
            cy + float(rotated[:, 1].min()) - margin,
            cy + float(rotated[:, 1].max()) + margin,
            cz + float(rotated[:, 2].min()) - margin,
            cz + float(rotated[:, 2].max()) + margin,
        )

    def generate(self: Box, **_: KwargsGenerateType) -> CadShape:
        """Generate a box CAD shape (OCCT).  Requires the ``[cad]`` extra."""
        from microgen.cad import make_box

        shape = make_box(self.dim, self.center)
        return rotate(shape, self.center, self.orientation)

    def generate_surface_mesh(
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
