from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
from typing import List, Tuple
import math as m
from scipy.spatial import KDTree

_UNIT_CUBE_SIZE = 1.0
_STRUT_NUMBER = 12
_STRUT_HEIGHTS = m.sqrt(2.0) / 2.0

_BASE_VERTICES = np.array([
            [sign * _UNIT_CUBE_SIZE/2.0 if i == axis else 0.0 for i in range(3)]
            for axis in range(3) for sign in [-1, 1]
        ])

def _compute_vertex_pairs() -> List[Tuple[int, int]]:
    tree = KDTree(_BASE_VERTICES)
    pairs = set()
    for i in range(len(_BASE_VERTICES)):
        _, indices = tree.query(_BASE_VERTICES[i], k=5)
        for j in indices[1:]:
            pairs.add(tuple(sorted((i, j))))
    return np.array(list(pairs))

_VERTEX_PAIRS = _compute_vertex_pairs()

class Octahedron(AbstractLattice):
    """
    Class to create a unit octahedron lattice of given cell size and density or strut radius
    """

    def __init__(self,
                 *args, **kwargs
                 ) -> None:
        super().__init__(*args, **kwargs, strut_number=_STRUT_NUMBER, strut_heights=_STRUT_HEIGHTS)
        
    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        """Compute the vertices of the octahedric lattice."""
        return self.center + self.cell_size * _BASE_VERTICES
    
    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        """Compute the centers of the struts."""
        return np.mean(self.vertices[_VERTEX_PAIRS], axis=1)
    
    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        """Compute the normalized direction vectors of the struts."""
        vectors = np.diff(self.vertices[_VERTEX_PAIRS], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)