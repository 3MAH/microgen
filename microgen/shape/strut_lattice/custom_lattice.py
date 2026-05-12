"""
===================================================================
Custom Lattice (:mod:`microgen.shape.strut_lattice.custom_lattice`)
===================================================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .abstract_lattice import AbstractLattice

if TYPE_CHECKING:
    import numpy as np
    import numpy.typing as npt
    from scipy.spatial.transform import Rotation

    from microgen.shape import Vector3DType


class CustomLattice(AbstractLattice):
    """Strut-based lattice from user-defined vertices and connectivity.

    ``base_vertices`` and ``strut_vertex_pairs`` are required positional.
    Configure the rest via chained ``with_*`` setters before calling
    :meth:`generate_cad` or :meth:`generate_surface_mesh`.
    """

    def __init__(
        self,
        base_vertices: npt.NDArray[np.float64],
        strut_vertex_pairs: npt.NDArray[np.int64],
        *,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType | Rotation = (0, 0, 0),
    ) -> None:
        """Initialize with the user-defined vertex layout and connectivity."""
        super().__init__(center=center, orientation=orientation)
        self._user_base_vertices = base_vertices
        self._user_strut_vertex_pairs = strut_vertex_pairs

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        # The base class accessor returns _user_base_vertices when set, so
        # this hook should never run for CustomLattice.
        return self._user_base_vertices

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        return self._user_strut_vertex_pairs
