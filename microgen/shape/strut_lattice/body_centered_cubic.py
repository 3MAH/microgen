import numpy as np
import numpy.typing as npt
import math as m
from itertools import product
from .abstract_lattice import AbstractLattice

class BodyCenteredCubic(AbstractLattice):
    """
    Class to create a unit body-centered cubic lattice of given cell size and density or strut radius
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, strut_heights=m.sqrt(3.0) / 2.0)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        unit_cube_vertices = np.array(list(product([-self._UNIT_CUBE_SIZE/2, self._UNIT_CUBE_SIZE/2], repeat=3)))        
        return np.vstack(([0.0, 0.0, 0.0], unit_cube_vertices))
    
    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        return np.array([[0, i] for i in range(1, 9)])