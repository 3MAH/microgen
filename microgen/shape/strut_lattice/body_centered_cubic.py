from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt
import math as m


class BodyCenteredCubic(AbstractLattice):
    """
    Class to create a unit body-centered cubic lattice of given cell size and density or strut radius
    """

    def __init__(self,
                 *args, **kwargs
                 ) -> None:

        super().__init__(*args, **kwargs, strut_number=8, strut_heights=m.sqrt(3.0)/2.0)
    
    ##TODO: reduce total number of struts

    def _compute_vertices(self) -> npt.NDArray[np.float64]:
        vertices_array = self.center + self.cell_size * np.array([
            [0.0, 0.0, 0.0],
            [0.5, -0.5, 0.5],
            [0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
            [-0.5, -0.5, 0.5],
            [0.5, -0.5, -0.5],
            [0.5, 0.5, -0.5],
            [-0.5, 0.5, -0.5],
            [-0.5, -0.5, -0.5]
        ])

        return vertices_array

    def _compute_strut_centers(self) -> npt.NDArray[np.float64]:
        centers_array = np.array([
            (self.vertices[1] + self.vertices[0]),
            (self.vertices[2] + self.vertices[0]),
            (self.vertices[3] + self.vertices[0]),
            (self.vertices[4] + self.vertices[0]),
            (self.vertices[5] + self.vertices[0]),
            (self.vertices[6] + self.vertices[0]),
            (self.vertices[7] + self.vertices[0]),
            (self.vertices[8] + self.vertices[0]),
        ]) / 2.0
        return centers_array

    def _compute_strut_directions(self) -> npt.NDArray[np.float64]:
        directions_array = np.array([
            (self.vertices[1] - self.vertices[0]) / np.linalg.norm((self.vertices[1] - self.vertices[0])),
            (self.vertices[2] - self.vertices[0]) / np.linalg.norm((self.vertices[2] - self.vertices[0])),
            (self.vertices[3] - self.vertices[0]) / np.linalg.norm((self.vertices[3] - self.vertices[0])),
            (self.vertices[4] - self.vertices[0]) / np.linalg.norm((self.vertices[4] - self.vertices[0])),
            (self.vertices[5] - self.vertices[0]) / np.linalg.norm((self.vertices[5] - self.vertices[0])),
            (self.vertices[6] - self.vertices[0]) / np.linalg.norm((self.vertices[6] - self.vertices[0])),
            (self.vertices[7] - self.vertices[0]) / np.linalg.norm((self.vertices[7] - self.vertices[0])),
            (self.vertices[8] - self.vertices[0]) / np.linalg.norm((self.vertices[8] - self.vertices[0])),
        ])

        return directions_array