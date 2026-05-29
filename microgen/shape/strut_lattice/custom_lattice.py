"""
===================================================================
Custom Lattice (:mod:`microgen.shape.strut_lattice.custom_lattice`)
===================================================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .abstract_lattice import AbstractLattice

if TYPE_CHECKING:
    import numpy.typing as npt
    from scipy.spatial.transform import Rotation

    from microgen.shape import Vector3DType

# Default topology: three orthogonal struts crossing at the unit-cube center.
_DEFAULT_BASE_VERTICES = np.array(
    [
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, -0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, -0.5],
        [0.0, 0.0, 0.5],
    ],
)
_DEFAULT_STRUT_VERTEX_PAIRS = np.array([[0, 1], [2, 3], [4, 5]])


class CustomLattice(AbstractLattice):
    """Strut-based lattice from user-defined vertices and connectivity.

    All parameters have defaults: the default topology is three orthogonal
    struts crossing at the unit-cube center (axis crosshair).  Configure
    everything via chained ``with_*`` setters before calling
    :meth:`generate_cad` or :meth:`generate_surface_mesh`.
    """

    _DEFAULT_STRUT_HEIGHTS = AbstractLattice._UNIT_CUBE_SIZE

    def __init__(
        self,
        base_vertices: npt.NDArray[np.float64] | None = None,
        strut_vertex_pairs: npt.NDArray[np.int64] | None = None,
        *,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType | Rotation = (0, 0, 0),
    ) -> None:
        """Initialize with the user-defined vertex layout and connectivity."""
        super().__init__(center=center, orientation=orientation)
        self._user_base_vertices = (
            base_vertices if base_vertices is not None else _DEFAULT_BASE_VERTICES
        )
        self._user_strut_vertex_pairs = (
            strut_vertex_pairs
            if strut_vertex_pairs is not None
            else _DEFAULT_STRUT_VERTEX_PAIRS
        )

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        # The base class accessor returns _user_base_vertices when set, so
        # this hook should never run for CustomLattice.
        return self._user_base_vertices

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        return self._user_strut_vertex_pairs
