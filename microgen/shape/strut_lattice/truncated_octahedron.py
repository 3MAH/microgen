from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
import math as m


class TruncatedOctahedron(AbstractLattice):
    """
    Class to create a unit truncated octahedron lattice of given cell size and density or strut radius
    """

    def __init__(self,
                 *args, **kwargs
                 ) -> None:

        super().__init__(*args, **kwargs, strut_number=36, strut_heights=m.sqrt(2.0) / 4.0)


    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        vertices_array = self.center + self.cell_size * np.array([
            [0.0, 1.0, 2.0],
            [0.0, -1.0, 2.0],
            [0.0, 1.0, -2.0],
            [0.0, -1.0, -2.0],
            [0.0, 2.0, 1.0],
            [0.0, -2.0, 1.0],
            [0.0, 2.0, -1.0],
            [0.0, -2.0, -1.0],
            [1.0, 0.0, 2.0],
            [1.0, 0.0, -2.0],
            [1.0, 2.0, 0.0],
            [1.0, -2.0, 0.0],
            [-1.0, 0.0, 2.0],
            [-1.0, 0.0, -2.0],
            [-1.0, 2.0, 0.0],
            [-1.0, -2.0, 0.0],
            [2.0, 0.0, 1.0],
            [2.0, 0.0, -1.0],
            [2.0, 1.0, 0.0],
            [2.0, -1.0, 0.0],
            [-2.0, 0.0, 1.0],
            [-2.0, 0.0, -1.0],
            [-2.0, 1.0, 0.0],
            [-2.0, -1.0, 0.0],
        ]) / 4.0

        return vertices_array

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        centers_array = np.array([
            (self.vertices[12] + self.vertices[0]),
            (self.vertices[1] + self.vertices[12]),
            (self.vertices[1] + self.vertices[8]),
            (self.vertices[8] + self.vertices[0]),
            (self.vertices[16] + self.vertices[8]),
            (self.vertices[16] + self.vertices[18]),
            (self.vertices[18] + self.vertices[10]),
            (self.vertices[10] + self.vertices[4]),
            (self.vertices[4] + self.vertices[0]),
            (self.vertices[5] + self.vertices[1]),
            (self.vertices[11] + self.vertices[5]),
            (self.vertices[19] + self.vertices[11]),
            (self.vertices[19] + self.vertices[16]),
            (self.vertices[17] + self.vertices[18]),
            (self.vertices[17] + self.vertices[19]),
            (self.vertices[6] + self.vertices[10]),
            (self.vertices[6] + self.vertices[2]),
            (self.vertices[9] + self.vertices[2]),
            (self.vertices[9] + self.vertices[17]),
            (self.vertices[3] + self.vertices[9]),
            (self.vertices[7] + self.vertices[3]),
            (self.vertices[11] + self.vertices[7]),
            (self.vertices[13] + self.vertices[2]),
            (self.vertices[13] + self.vertices[3]),
            (self.vertices[21] + self.vertices[13]),
            (self.vertices[23] + self.vertices[21]),
            (self.vertices[15] + self.vertices[23]),
            (self.vertices[7] + self.vertices[15]),
            (self.vertices[5] + self.vertices[15]),
            (self.vertices[20] + self.vertices[23]),
            (self.vertices[12] + self.vertices[20]),
            (self.vertices[22] + self.vertices[20]),
            (self.vertices[21] + self.vertices[22]),
            (self.vertices[14] + self.vertices[22]),
            (self.vertices[14] + self.vertices[4]),
            (self.vertices[6] + self.vertices[14])
        ]) / 2.0
        return centers_array

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        directions_array = np.array([
            (self.vertices[12] - self.vertices[0]) / np.linalg.norm((self.vertices[12] - self.vertices[0])),
            (self.vertices[1] - self.vertices[12]) / np.linalg.norm((self.vertices[1] - self.vertices[12])),
            (self.vertices[1] - self.vertices[8]) / np.linalg.norm((self.vertices[1] - self.vertices[8])),
            (self.vertices[8] - self.vertices[0]) / np.linalg.norm((self.vertices[8] - self.vertices[0])),
            (self.vertices[16] - self.vertices[8]) / np.linalg.norm((self.vertices[16] - self.vertices[8])),
            (self.vertices[16] - self.vertices[18]) / np.linalg.norm((self.vertices[16] - self.vertices[18])),
            (self.vertices[18] - self.vertices[10]) / np.linalg.norm((self.vertices[18] - self.vertices[10])),
            (self.vertices[10] - self.vertices[4]) / np.linalg.norm((self.vertices[10] - self.vertices[4])),
            (self.vertices[4] - self.vertices[0]) / np.linalg.norm((self.vertices[4] - self.vertices[0])),
            (self.vertices[5] - self.vertices[1]) / np.linalg.norm((self.vertices[5] - self.vertices[1])),
            (self.vertices[11] - self.vertices[5]) / np.linalg.norm((self.vertices[11] - self.vertices[5])),
            (self.vertices[19] - self.vertices[11]) / np.linalg.norm((self.vertices[19] - self.vertices[11])),
            (self.vertices[19] - self.vertices[16]) / np.linalg.norm((self.vertices[19] - self.vertices[16])),
            (self.vertices[17] - self.vertices[18]) / np.linalg.norm((self.vertices[17] - self.vertices[18])),
            (self.vertices[17] - self.vertices[19]) / np.linalg.norm((self.vertices[17] - self.vertices[19])),
            (self.vertices[6] - self.vertices[10]) / np.linalg.norm((self.vertices[6] - self.vertices[10])),
            (self.vertices[6] - self.vertices[2]) / np.linalg.norm((self.vertices[6] - self.vertices[2])),
            (self.vertices[9] - self.vertices[2]) / np.linalg.norm((self.vertices[9] - self.vertices[2])),
            (self.vertices[9] - self.vertices[17]) / np.linalg.norm((self.vertices[9] - self.vertices[17])),
            (self.vertices[3] - self.vertices[9]) / np.linalg.norm((self.vertices[3] - self.vertices[9])),
            (self.vertices[7] - self.vertices[3]) / np.linalg.norm((self.vertices[7] - self.vertices[3])),
            (self.vertices[11] - self.vertices[7]) / np.linalg.norm((self.vertices[11] - self.vertices[7])),
            (self.vertices[13] - self.vertices[2]) / np.linalg.norm((self.vertices[13] - self.vertices[2])),
            (self.vertices[13] - self.vertices[3]) / np.linalg.norm((self.vertices[13] - self.vertices[3])),
            (self.vertices[21] - self.vertices[13]) / np.linalg.norm((self.vertices[21] - self.vertices[13])),
            (self.vertices[23] - self.vertices[21]) / np.linalg.norm((self.vertices[23] - self.vertices[21])),
            (self.vertices[15] - self.vertices[23]) / np.linalg.norm((self.vertices[15] - self.vertices[23])),
            (self.vertices[7] - self.vertices[15]) / np.linalg.norm((self.vertices[7] - self.vertices[15])),
            (self.vertices[5] - self.vertices[15]) / np.linalg.norm((self.vertices[5] - self.vertices[15])),
            (self.vertices[20] - self.vertices[23]) / np.linalg.norm((self.vertices[20] - self.vertices[23])),
            (self.vertices[12] - self.vertices[20]) / np.linalg.norm((self.vertices[12] - self.vertices[20])),
            (self.vertices[22] - self.vertices[20]) / np.linalg.norm((self.vertices[22] - self.vertices[20])),
            (self.vertices[21] - self.vertices[22]) / np.linalg.norm((self.vertices[21] - self.vertices[22])),
            (self.vertices[14] - self.vertices[22]) / np.linalg.norm((self.vertices[14] - self.vertices[22])),
            (self.vertices[14] - self.vertices[4]) / np.linalg.norm((self.vertices[14] - self.vertices[4])),
            (self.vertices[6] - self.vertices[14]) / np.linalg.norm((self.vertices[6] - self.vertices[14]))
        ])

        return directions_array