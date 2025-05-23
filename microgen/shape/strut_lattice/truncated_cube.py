"""
===================================================================
Truncated Cube (:mod:`microgen.shape.strut_lattice.truncated_cube`)
===================================================================
"""

from itertools import product

import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree

from .abstract_lattice import AbstractLattice


class TruncatedCube(AbstractLattice):
    """
    Class to create a unit truncated cubic lattice of given cell size and density or strut radius

    .. jupyter-execute::
       :hide-code:
       :hide-output:

       import microgen

       shape = microgen.TruncatedCube(strut_radius=0.1).generate_vtk()

    .. jupyter-execute::
       :hide-code:

       shape.plot(color='white')
    """

    def __init__(self, *args, **kwargs) -> None:
        kwargs.setdefault("strut_heights", np.sqrt(2.0) - 1.0)
        super().__init__(*args, **kwargs)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        values = [
            -self._UNIT_CUBE_SIZE / 2.0 * (np.sqrt(2) - 1.0),
            self._UNIT_CUBE_SIZE / 2.0 * (np.sqrt(2) - 1.0),
            -self._UNIT_CUBE_SIZE / 2.0,
            self._UNIT_CUBE_SIZE / 2.0,
        ]

        vertices = set()
        for x, y, z in product(values, repeat=3):
            if (
                sum(abs(coord) == self._UNIT_CUBE_SIZE / 2.0 for coord in (x, y, z))
                == 2
            ):
                vertices.add((x, y, z))

        return np.array(list(vertices))

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        tree = KDTree(self.base_vertices)
        pairs = set()
        TOLERANCE = 1e-5
        distance = np.sqrt(2) - 1.0 + TOLERANCE

        for i, vertex in enumerate(self.base_vertices):
            indices = tree.query_ball_point(vertex, distance)
            for j in indices:
                if i < j:
                    pairs.add((i, j))

        return np.array(list(pairs))
