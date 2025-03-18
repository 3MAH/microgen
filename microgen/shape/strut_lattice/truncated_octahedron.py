from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
import math as m
from itertools import product, permutations
from scipy.spatial import KDTree

class TruncatedOctahedron(AbstractLattice):
    """
    Class to create a unit truncated octahedron lattice of given cell size and density or strut radius
    """

    def __init__(self,
                 *args, **kwargs
                 ) -> None:
        self._base_vertices = self._generate_base_vertices()
        self._strut_vertex_pairs = self._generate_strut_vertex_pairs()
        super().__init__(*args, **kwargs, strut_number=36, strut_heights=m.sqrt(2.0) / 4.0)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        base_vertices = set()
        for perm in permutations([0, self._UNIT_CUBE_SIZE/4.0, self._UNIT_CUBE_SIZE/2.0]):
            for signs in product([-1, 1], repeat=3):
                vertex = tuple(s * p for s, p in zip(signs, perm))
                if sum(abs(v) for v in vertex) == 3.0/4.0:
                    base_vertices.add(vertex)
        return np.array(list(base_vertices))

    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        return self.center + self.cell_size * self._base_vertices
    
    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        """Generate index pairs representing the struts in the truncated octahedron."""
        tree = KDTree(self._base_vertices)
        tolerance = 1e-5
        pairs = set()
        connection_distance = m.sqrt(2.0) / 4.0 + tolerance
        for i, vertex in enumerate(self._base_vertices):
            indices = tree.query_ball_point(vertex, connection_distance)
            for j in indices:
                if i < j:
                    pairs.add((i, j))
        return np.array(list(pairs))

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        return np.mean(self.vertices[self._strut_vertex_pairs], axis=1)

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        vectors = np.diff(self.vertices[self._strut_vertex_pairs], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)