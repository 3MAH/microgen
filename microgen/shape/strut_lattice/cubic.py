import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree
from itertools import product
from .abstract_lattice import AbstractLattice


class Cubic(AbstractLattice):
    """
    Class to create a unit cubic lattice of given cell size and density or strut radius.
    """

    def __init__(self, *args, **kwargs) -> None:
        self._base_vertices = self._generate_base_vertices()
        self._strut_vertex_pairs = self._generate_vertex_pairs()
        super().__init__(*args, **kwargs, strut_number=12, strut_heights=self._UNIT_CUBE_SIZE)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        return np.array(list(product([-self._UNIT_CUBE_SIZE/2, self._UNIT_CUBE_SIZE/2], repeat=3)))

    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        return self.center + self.cell_size * self._base_vertices
    
    def _generate_vertex_pairs(self) -> npt.NDArray[np.int64]:
        return np.array([
            [i, j] for i in range(len(self._base_vertices)) 
            for j in KDTree(self._base_vertices).query_ball_point(self._base_vertices[i], r=self._UNIT_CUBE_SIZE) if i < j
        ])

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        return np.mean(self.vertices[self._strut_vertex_pairs], axis=1)

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        vectors = np.diff(self.vertices[self._strut_vertex_pairs], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)

