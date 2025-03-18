from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
import math as m
from itertools import product
from scipy.spatial import KDTree

class RhombicDodecahedron(AbstractLattice):
    """
    Class to create a unit rhombic dodecahedron lattice of given cell size and density or strut radius
    """

    def __init__(self,
                 *args, **kwargs
                 ) -> None:
        self._base_vertices = self._generate_base_vertices()
        self._strut_vertex_pairs = self._generate_strut_vertex_pairs()
        super().__init__(*args, **kwargs, strut_number=32, strut_heights=m.sqrt(3.0) / 4.0)
        
    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        outer_cube_vertices = list(product([-self._UNIT_CUBE_SIZE/2, self._UNIT_CUBE_SIZE/2], repeat=3))
        inner_cube_vertices = list(product([-self._UNIT_CUBE_SIZE/4, self._UNIT_CUBE_SIZE/4], repeat=3))
                
        face_centers = [
            [sign * self._UNIT_CUBE_SIZE/2 if i == axis else 0.0 for i in range(3)]
            for axis in range(3) for sign in [-1, 1]
        ]
        
        return np.array(outer_cube_vertices + inner_cube_vertices + face_centers)
    
    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        return self.center + self.cell_size * self._base_vertices
    
    def _generate_strut_vertex_pairs(self) -> npt.NDArray[int]:
        tree = KDTree(self._base_vertices)
        tolerance = 1e-5
        target_distance = m.sqrt(3.0) / 4.0 + tolerance
        pairs = set(tree.query_pairs(r=target_distance))
        return np.array(list(pairs))
    
    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        return np.mean(self.vertices[self._strut_vertex_pairs], axis=1)
    
    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        vectors = np.diff(self.vertices[self._strut_vertex_pairs], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
        