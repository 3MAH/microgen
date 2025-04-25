"""
=================================================
Cubic (:mod:`microgen.shape.strut_lattice.cubic`)
=================================================
"""

from itertools import product

import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree

from .abstract_lattice import AbstractLattice


class Cubic(AbstractLattice):
    """
    Class to create a unit cubic lattice of given cell size and density or strut radius.

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Cubic(strut_radius=0.1).generate_vtk()
       shape.plot(color='white')
    """

    def __init__(self, *args, **kwargs) -> None:
        kwargs.setdefault("strut_heights", self._UNIT_CUBE_SIZE)
        super().__init__(*args, **kwargs)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        return np.array(
            list(
                product([-self._UNIT_CUBE_SIZE / 2, self._UNIT_CUBE_SIZE / 2], repeat=3)
            )
        )

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        return np.array(
            [
                [i, j]
                for i in range(len(self.base_vertices))
                for j in KDTree(self.base_vertices).query_ball_point(
                    self.base_vertices[i], r=self._UNIT_CUBE_SIZE
                )
                if i < j
            ]
        )
