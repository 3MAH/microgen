from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
import math as m
from scipy.spatial import KDTree

class Octahedron(AbstractLattice):
    """
    Class to create a unit octahedron lattice of given cell size and density or strut radius
    """

    def __init__(self,
                 *args, **kwargs
                 ) -> None:
        self._base_vertices = self._generate_base_vertices()
        self._strut_vertex_pairs = self._generate_vertex_pairs()
        super().__init__(*args, **kwargs, strut_number=12, strut_heights=m.sqrt(2.0) / 2.0)
    
    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        return np.array([
            [sign * self._UNIT_CUBE_SIZE/2.0 if i == axis else 0.0 for i in range(3)]
            for axis in range(3) for sign in [-1, 1]
        ])
        
    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        """Compute the vertices of the octahedric lattice."""
        return self.center + self.cell_size * self._base_vertices
    
    def _generate_vertex_pairs(self) -> npt.NDArray[int]:
        tree = KDTree(self._base_vertices)
        pairs = set()
        for i in range(len(self._base_vertices)):
            _, indices = tree.query(self._base_vertices[i], k=5)
            for j in indices[1:]:
                pairs.add(tuple(sorted((i, j))))
        return np.array(list(pairs))
    
    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        """Compute the centers of the struts."""
        return np.mean(self.vertices[self._strut_vertex_pairs], axis=1)
    
    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        """Compute the normalized direction vectors of the struts."""
        vectors = np.diff(self.vertices[self._strut_vertex_pairs], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)