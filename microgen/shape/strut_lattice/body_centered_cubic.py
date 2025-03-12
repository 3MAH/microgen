import numpy as np
import numpy.typing as npt
import math as m
from itertools import product
from .abstract_lattice import AbstractLattice

_UNIT_CUBE_SIZE = 1.0
_STRUT_NUMBER = 8
_STRUT_HEIGHTS = m.sqrt(3.0) / 2.0

class BodyCenteredCubic(AbstractLattice):
    """
    Class to create a unit body-centered cubic lattice of given cell size and density or strut radius
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, strut_number=_STRUT_NUMBER, strut_heights=_STRUT_HEIGHTS)

    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        """Compute the vertices of the BCC lattice."""
        unit_cube_vertices = np.array(list(product([-_UNIT_CUBE_SIZE/2, _UNIT_CUBE_SIZE/2], repeat=3)))
        vertices = np.vstack(([0.0, 0.0, 0.0], unit_cube_vertices))
        return self.center + self.cell_size * vertices

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        """Compute the centers of the struts."""
        return (self.vertices[1:] + self.vertices[0]) / 2.0

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        """Compute the normalized direction vectors of the struts."""
        vectors = self.vertices[1:] - self.vertices[0]
        return vectors / np.linalg.norm(vectors, axis=1)[:, np.newaxis]
