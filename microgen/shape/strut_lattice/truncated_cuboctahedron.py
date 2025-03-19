from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
import math as m
from itertools import product, permutations
from scipy.spatial import KDTree


class TruncatedCuboctahedron(AbstractLattice):
    """
    Class to create a unit truncated cuboctahedron lattice of given cell size and density or strut radius
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, strut_heights=1 / (1.0 + 2.0 * m.sqrt(2)))

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        values = [1, -(1 + m.sqrt(2)), -(1 + 2 * m.sqrt(2))]
        signs = [1, -1]
        result = []

        for perm in permutations(values):
            for sign in product(signs, repeat=3):
                result.append([sign[0] * perm[0], sign[1] * perm[1], sign[2] * perm[2]])

        return np.array(result) / (2.0 + 4.0 * m.sqrt(2))

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        kdtree = KDTree(self.base_vertices)
        tolerance = 1e-5
        threshold_distance = 1 / (1.0 + 2.0 * m.sqrt(2)) + tolerance

        pairs = []
        for i in range(len(self.base_vertices)):
            neighbors = kdtree.query_ball_point(
                self.base_vertices[i], threshold_distance
            )
            for j in neighbors:
                if i < j:
                    pairs.append((i, j))

        return np.array(pairs)
