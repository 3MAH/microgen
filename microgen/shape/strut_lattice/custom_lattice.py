"""
===================================================================
Custom Lattice (:mod:`microgen.shape.strut_lattice.custom_lattice`)
===================================================================
"""

import numpy as np
import numpy.typing as npt

from .abstract_lattice import AbstractLattice


class CustomLattice(AbstractLattice):
    """
    Class to create a custom lattice with user-defined base vertices and strut vertex pairs
    """

    def __init__(
        self,
        base_vertices: npt.NDArray[np.float64],
        strut_vertex_pairs: npt.NDArray[np.int64],
        *args,
        **kwargs,
    ) -> None:
        kwargs["base_vertices"] = base_vertices
        kwargs["strut_vertex_pairs"] = strut_vertex_pairs
        super().__init__(
            *args,
            **kwargs,
        )

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        return self.base_vertices

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        return self.strut_vertex_pairs
