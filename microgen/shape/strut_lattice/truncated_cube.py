from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
import math as m
from itertools import product
from scipy.spatial import KDTree


class TruncatedCube(AbstractLattice):
    """
    Class to create a unit truncated cubic lattice of given cell size and density or strut radius
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, strut_heights=m.sqrt(2) - 1.0)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        values = [
            -self._UNIT_CUBE_SIZE / 2.0 * (m.sqrt(2) - 1.0),
            self._UNIT_CUBE_SIZE / 2.0 * (m.sqrt(2) - 1.0),
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
        tolerance = 1e-5
        distance = m.sqrt(2) - 1.0 + tolerance

        for i, vertex in enumerate(self.base_vertices):
            indices = tree.query_ball_point(vertex, distance)
            for j in indices:
                if i < j:
                    pairs.add((i, j))

        return np.array(list(pairs))
