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


class CustomLattice(AbstractLattice):
    """Class to create a custom lattice.

    Uses user-defined base vertices and strut vertex pairs.
    """

    def __init__(  # noqa: PLR0913
        self,
        base_vertices: npt.NDArray[np.float64],
        strut_vertex_pairs: npt.NDArray[np.int64],
        strut_radius: float | None = None,
        strut_heights: float | list[float] | None = None,
        cell_size: float = 1.0,
        strut_joints: bool = False,
        density: float | None = None,
    ) -> None:
        """Initialize the custom lattice."""
        super().__init__(
            base_vertices=base_vertices,
            strut_vertex_pairs=strut_vertex_pairs,
            strut_radius=strut_radius,
            strut_heights=strut_heights,
            cell_size=cell_size,
            strut_joints=strut_joints,
            density=density,
        )

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        return self.base_vertices

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        return self.strut_vertex_pairs
