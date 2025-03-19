from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
import math as m
from itertools import product, permutations
from scipy.spatial import KDTree


class RhombicCuboctahedron(AbstractLattice):
    """
    Class to create a unit rhombic cuboctahedron lattice of given cell size and density or strut radius
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, strut_heights=m.sqrt(2.0) - 1.0)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        permutations_set = set(
            permutations(
                [self._UNIT_CUBE_SIZE / 2.0, (np.sqrt(2) - 1) / 2, (np.sqrt(2) - 1) / 2]
            )
        )

        vertices = []
        for permutation in permutations_set:
            for signs in product([-1, 1], repeat=3):
                vertex = tuple(s * p for s, p in zip(signs, permutation))
                vertices.append(vertex)
        return np.array(vertices)

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
