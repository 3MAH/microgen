from .abstract_lattice import AbstractLattice
from itertools import product
import numpy as np
import numpy.typing as npt
import math as m
from scipy.spatial import KDTree

class OctetTruss(AbstractLattice):
    """
    Class to create a unit octet-truss lattice of given cell size and density or strut radius
    """

    def __init__(self,
                 *args, **kwargs
                 ) -> None:
        self._base_vertices = self._generate_base_vertices()
        self._strut_vertex_pairs = self._generate_strut_vertex_pairs()
        super().__init__(*args, **kwargs, strut_number=36, strut_heights=m.sqrt(2.0) / 2.0)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        cube_vertices = list(product([-self._UNIT_CUBE_SIZE/2, self._UNIT_CUBE_SIZE/2], repeat=3))
        
        face_centers = [
            [sign * self._UNIT_CUBE_SIZE/2 if i == axis else 0.0 for i in range(3)]
            for axis in range(3) for sign in [-1, 1]
        ]

        return np.array(cube_vertices + face_centers)
    
    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        """Compute the vertices of the face-centered cubic lattice."""
        return self.center + self.cell_size * self._base_vertices
    
    def _generate_strut_vertex_pairs(self) -> npt.NDArray[int]:
        """Generate index pairs representing the struts in the octet-truss using KDTree."""
        tree = KDTree(self._base_vertices)
        pairs = set()
        tolerance = 1e-5
        
        # Define a reasonable connection distance threshold
        connection_distance = (self._UNIT_CUBE_SIZE / m.sqrt(2)) + tolerance
        
        for i, vertex in enumerate(self._base_vertices):
            indices = tree.query_ball_point(vertex, connection_distance)
            for j in indices:
                if i != j:
                    pairs.add(tuple(sorted((i, j))))
        
        return np.array(list(pairs))

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        return np.mean(self.vertices[self._strut_vertex_pairs], axis=1)

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        vectors = np.diff(self.vertices[self._strut_vertex_pairs], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)