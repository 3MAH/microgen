import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree
import itertools
from .abstract_lattice import AbstractLattice

_UNIT_CUBE_SIZE = 1.0
_STRUT_NUMBER = 12
_STRUT_HEIGHTS = 1.0
_VERTICES = np.array(list(itertools.product([-_UNIT_CUBE_SIZE/2, _UNIT_CUBE_SIZE/2], repeat=3)))
_STRUT_VERTEX_PAIRS = np.array([
    [i, j] for i in range(len(_VERTICES)) 
    for j in KDTree(_VERTICES).query_ball_point(_VERTICES[i], r=_STRUT_HEIGHTS) if i < j
])

class Cubic(AbstractLattice):
    """
    Class to create a unit cubic lattice of given cell size and density or strut radius.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, strut_number=_STRUT_NUMBER, strut_heights=_STRUT_HEIGHTS)

    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        """Compute the vertices of the cubic lattice."""
        return self.center + self.cell_size * _VERTICES

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        """Compute the centers of the struts."""
        return np.mean(self.vertices[_STRUT_VERTEX_PAIRS], axis=1)

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        """Compute the normalized direction vectors of the struts."""
        vectors = np.diff(self.vertices[_STRUT_VERTEX_PAIRS], axis=1).squeeze()
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)
