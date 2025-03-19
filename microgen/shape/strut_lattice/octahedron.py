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
        super().__init__(*args, **kwargs, strut_heights=m.sqrt(2.0) / 2.0)
    
    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        return np.array([
            [sign * self._UNIT_CUBE_SIZE/2.0 if i == axis else 0.0 for i in range(3)]
            for axis in range(3) for sign in [-1, 1]
        ])
    
    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        tree = KDTree(self.base_vertices)
        pairs = set()
        for i in range(len(self.base_vertices)):
            _, indices = tree.query(self.base_vertices[i], k=5)
            for j in indices[1:]:
                pairs.add(tuple(sorted((i, j))))
        return np.array(list(pairs))