from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
import math as m

_STRUT_NUMBER = 12
_STRUT_HEIGHTS = m.sqrt(2.0) / 2.0

class Octahedron(AbstractLattice):
    """
    Class to create a unit octahedron lattice of given cell size and density or strut radius
    """

    def __init__(self,
                 *args, **kwargs
                 ) -> None:
        super().__init__(*args, **kwargs, strut_number=_STRUT_NUMBER, strut_heights=_STRUT_HEIGHTS)
        
    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        vertices_array = self.center + self.cell_size * np.array([
            [0.5, 0.0, 0.0],
            [-0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, -0.5, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, -0.5]
        ])

        return vertices_array
    
    def _compute_strut_centers(self):
        centers_array = np.array([
            (self.vertices[0] + self.vertices[2]),
            (self.vertices[0] + self.vertices[3]),
            (self.vertices[0] + self.vertices[4]),
            (self.vertices[0] + self.vertices[5]),
            (self.vertices[1] + self.vertices[2]),
            (self.vertices[1] + self.vertices[3]),
            (self.vertices[1] + self.vertices[4]),
            (self.vertices[1] + self.vertices[5]),
            (self.vertices[2] + self.vertices[4]),
            (self.vertices[2] + self.vertices[5]),
            (self.vertices[3] + self.vertices[4]),
            (self.vertices[3] + self.vertices[5])
        ]) / 2.0
        return centers_array
    
    def _compute_strut_directions(self):
        directions_array = np.array([
            (self.vertices[2] - self.vertices[0]) / np.linalg.norm((self.vertices[2] - self.vertices[0])),
            (self.vertices[3] - self.vertices[0]) / np.linalg.norm((self.vertices[3] - self.vertices[0])),
            (self.vertices[4] - self.vertices[0]) / np.linalg.norm((self.vertices[4] - self.vertices[0])),
            (self.vertices[5] - self.vertices[0]) / np.linalg.norm((self.vertices[5] - self.vertices[0])),
            (self.vertices[2] - self.vertices[1]) / np.linalg.norm((self.vertices[2] - self.vertices[1])),
            (self.vertices[3] - self.vertices[1]) / np.linalg.norm((self.vertices[3] - self.vertices[1])),
            (self.vertices[4] - self.vertices[1]) / np.linalg.norm((self.vertices[4] - self.vertices[1])),
            (self.vertices[5] - self.vertices[1]) / np.linalg.norm((self.vertices[5] - self.vertices[1])),
            (self.vertices[4] - self.vertices[2]) / np.linalg.norm((self.vertices[4] - self.vertices[2])),
            (self.vertices[5] - self.vertices[2]) / np.linalg.norm((self.vertices[5] - self.vertices[2])),
            (self.vertices[4] - self.vertices[3]) / np.linalg.norm((self.vertices[4] - self.vertices[3])),
            (self.vertices[5] - self.vertices[3]) / np.linalg.norm((self.vertices[5] - self.vertices[3]))
        ])
        return directions_array