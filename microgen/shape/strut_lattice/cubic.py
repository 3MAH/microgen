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
        self.base_vertices = np.array(list(product([-self._UNIT_CUBE_SIZE/2, self._UNIT_CUBE_SIZE/2], repeat=3)))
        self.strut_vertex_pairs = self._generate_vertex_pairs()
        super().__init__(*args, **kwargs, strut_number=12, strut_heights=self._UNIT_CUBE_SIZE)

    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        """Compute the vertices of the cubic lattice."""
        return self.center + self.cell_size * self.base_vertices

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        """Compute the centers of the struts."""
        return np.mean(self.vertices[self.strut_vertex_pairs], axis=1)

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        """Compute the normalized direction vectors of the struts."""
        vectors = np.diff(self.vertices[self.strut_vertex_pairs], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
    
    def _generate_vertex_pairs(self) -> npt.NDArray[int]:
        """
        Compute the strut vertex pairs.
        """
        return np.array([
            [i, j] for i in range(len(self.base_vertices)) 
            for j in KDTree(self.base_vertices).query_ball_point(self.base_vertices[i], r=self._UNIT_CUBE_SIZE) if i < j
        ])
