from .abstract_lattice import AbstractLattice
from itertools import product
import numpy as np
import numpy.typing as npt
import math as m
from scipy.spatial import KDTree

class Cuboctahedron(AbstractLattice):
    """
    Class to create a unit cuboctahedron lattice of given cell size and density or strut radius
    """

    def __init__(self,
                 *args, **kwargs
                 ) -> None:
        self._base_vertices = self._generate_base_vertices()
        self._strut_vertex_pairs = self._generate_vertex_pairs()
        super().__init__(*args, **kwargs, strut_number=24, strut_heights=m.sqrt(2.0)/2.0)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        cube_vertices = np.array(list(product([-self._UNIT_CUBE_SIZE/2, self._UNIT_CUBE_SIZE/2], repeat=3)))
        
        edges = [(i, j) for i in range(len(cube_vertices)) for j in range(i + 1, len(cube_vertices))
                 if np.sum(np.abs(cube_vertices[i] - cube_vertices[j])) == self._UNIT_CUBE_SIZE]
        
        return np.array([(cube_vertices[i] + cube_vertices[j]) / 2.0 for i, j in edges])

    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        vertices_array = self.center + self.cell_size * self._base_vertices

        return vertices_array
    
    def _generate_vertex_pairs(self) -> npt.NDArray[int]:
        tree = KDTree(self._base_vertices)
        pairs = set()
        tolerance = 1e-5
        for i, indices in enumerate(tree.query_ball_point(self._base_vertices, r=self._UNIT_CUBE_SIZE/m.sqrt(2.0) + tolerance)):
            for j in indices:
                if i != j:
                    pairs.add(tuple(sorted((i, j))))
        
        return np.array(list(pairs))

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        return np.mean(self.vertices[self._strut_vertex_pairs], axis=1)

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        vectors = np.diff(self.vertices[self._strut_vertex_pairs], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)