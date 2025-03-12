from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
import math as m

_STRUT_NUMBER = 32
_STRUT_HEIGHTS = m.sqrt(3.0) / 4.0

class RhombicDodecahedron(AbstractLattice):
    """
    Class to create a unit rhombic dodecahedron lattice of given cell size and density or strut radius
    """

    def __init__(self,
                 *args, **kwargs
                 ) -> None:
        super().__init__(*args, **kwargs, strut_number=_STRUT_NUMBER, strut_heights=_STRUT_HEIGHTS)
        
    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        vertices_array = self.center + self.cell_size * np.array([
            [0.25, 0.25, 0.25],
            [0.25, 0.25, -0.25],
            [0.25, -0.25, 0.25],
            [0.25, -0.25, -0.25],
            [-0.25, 0.25, 0.25],
            [-0.25, 0.25, -0.25],
            [-0.25, -0.25, 0.25],
            [-0.25, -0.25, -0.25],
            [0.5, 0.0, 0.0],
            [-0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, -0.5, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, -0.5],
            [0.5, 0.5, 0.5],
            [0.5, 0.5, -0.5],
            [0.5, -0.5, 0.5],
            [0.5, -0.5, -0.5],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, 0.5],
            [-0.5, -0.5, -0.5]
        ])

        return vertices_array
    
    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        centers_array = np.array([
            (self.vertices[0] + self.vertices[8]),
            (self.vertices[0] + self.vertices[10]),
            (self.vertices[0] + self.vertices[12]),
            (self.vertices[0] + self.vertices[14]),
            (self.vertices[1] + self.vertices[8]),
            (self.vertices[1] + self.vertices[10]),
            (self.vertices[1] + self.vertices[13]),
            (self.vertices[1] + self.vertices[15]),
            (self.vertices[2] + self.vertices[8]),
            (self.vertices[2] + self.vertices[11]),
            (self.vertices[2] + self.vertices[12]),
            (self.vertices[2] + self.vertices[16]),
            (self.vertices[3] + self.vertices[8]),
            (self.vertices[3] + self.vertices[11]),
            (self.vertices[3] + self.vertices[13]),
            (self.vertices[3] + self.vertices[17]),
            (self.vertices[4] + self.vertices[9]),
            (self.vertices[4] + self.vertices[10]),
            (self.vertices[4] + self.vertices[12]),
            (self.vertices[4] + self.vertices[18]),
            (self.vertices[5] + self.vertices[9]),
            (self.vertices[5] + self.vertices[10]),
            (self.vertices[5] + self.vertices[13]),
            (self.vertices[5] + self.vertices[19]),
            (self.vertices[6] + self.vertices[9]),
            (self.vertices[6] + self.vertices[11]),
            (self.vertices[6] + self.vertices[12]),
            (self.vertices[6] + self.vertices[20]),
            (self.vertices[7] + self.vertices[9]),
            (self.vertices[7] + self.vertices[11]),
            (self.vertices[7] + self.vertices[13]),
            (self.vertices[7] + self.vertices[21])
        ])/2.0
        return centers_array
    
    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        directions_array = np.array([
            (self.vertices[8] - self.vertices[0]) / np.linalg.norm((self.vertices[8] - self.vertices[0])),
            (self.vertices[10] - self.vertices[0]) / np.linalg.norm((self.vertices[10] - self.vertices[0])),
            (self.vertices[12] - self.vertices[0]) / np.linalg.norm((self.vertices[12] - self.vertices[0])),
            (self.vertices[14] - self.vertices[0]) / np.linalg.norm((self.vertices[14] - self.vertices[0])),
            (self.vertices[8] - self.vertices[1]) / np.linalg.norm((self.vertices[8] - self.vertices[1])),
            (self.vertices[10] - self.vertices[1]) / np.linalg.norm((self.vertices[10] - self.vertices[1])),
            (self.vertices[13] - self.vertices[1]) / np.linalg.norm((self.vertices[13] - self.vertices[1])),
            (self.vertices[15] - self.vertices[1]) / np.linalg.norm((self.vertices[15] - self.vertices[1])),
            (self.vertices[8] - self.vertices[2]) / np.linalg.norm((self.vertices[8] - self.vertices[2])),
            (self.vertices[11] - self.vertices[2]) / np.linalg.norm((self.vertices[11] - self.vertices[2])),
            (self.vertices[12] - self.vertices[2]) / np.linalg.norm((self.vertices[12] - self.vertices[2])),
            (self.vertices[16] - self.vertices[2]) / np.linalg.norm((self.vertices[16] - self.vertices[2])),
            (self.vertices[8] - self.vertices[3]) / np.linalg.norm((self.vertices[8] - self.vertices[3])),
            (self.vertices[11] - self.vertices[3]) / np.linalg.norm((self.vertices[11] - self.vertices[3])),
            (self.vertices[13] - self.vertices[3]) / np.linalg.norm((self.vertices[13] - self.vertices[3])),
            (self.vertices[17] - self.vertices[3]) / np.linalg.norm((self.vertices[17] - self.vertices[3])),
            (self.vertices[9] - self.vertices[4]) / np.linalg.norm((self.vertices[9] - self.vertices[4])),
            (self.vertices[10] - self.vertices[4]) / np.linalg.norm((self.vertices[10] - self.vertices[4])),
            (self.vertices[12] - self.vertices[4]) / np.linalg.norm((self.vertices[12] - self.vertices[4])),
            (self.vertices[18] - self.vertices[4]) / np.linalg.norm((self.vertices[18] - self.vertices[4])),
            (self.vertices[9] - self.vertices[5]) / np.linalg.norm((self.vertices[9] - self.vertices[5])),
            (self.vertices[10] - self.vertices[5]) / np.linalg.norm((self.vertices[10] - self.vertices[5])),
            (self.vertices[13] - self.vertices[5]) / np.linalg.norm((self.vertices[13] - self.vertices[5])),
            (self.vertices[19] - self.vertices[5]) / np.linalg.norm((self.vertices[19] - self.vertices[5])),
            (self.vertices[9] - self.vertices[6]) / np.linalg.norm((self.vertices[9] - self.vertices[6])),
            (self.vertices[11] - self.vertices[6]) / np.linalg.norm((self.vertices[11] - self.vertices[6])),
            (self.vertices[12] - self.vertices[6]) / np.linalg.norm((self.vertices[12] - self.vertices[6])),
            (self.vertices[20] - self.vertices[6]) / np.linalg.norm((self.vertices[20] - self.vertices[6])),
            (self.vertices[9] - self.vertices[7]) / np.linalg.norm((self.vertices[9] - self.vertices[7])),
            (self.vertices[11] - self.vertices[7]) / np.linalg.norm((self.vertices[11] - self.vertices[7])),
            (self.vertices[13] - self.vertices[7]) / np.linalg.norm((self.vertices[13] - self.vertices[7])),
            (self.vertices[21] - self.vertices[7]) / np.linalg.norm((self.vertices[21] - self.vertices[7]))
        ])
        return directions_array
        