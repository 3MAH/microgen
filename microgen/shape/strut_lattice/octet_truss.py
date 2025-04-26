"""
=============================================================
Octet-truss (:mod:`microgen.shape.strut_lattice.octet_truss`)
=============================================================
"""

from itertools import product

import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree

from .abstract_lattice import AbstractLattice


class OctetTruss(AbstractLattice):
    """
    Class to create a unit octet-truss lattice of given cell size and density or strut radius

    .. jupyter-execute::
       :hide-code:
       :hide-output:

       import microgen

       shape = microgen.OctetTruss(strut_radius=0.1).generate_vtk()

    .. jupyter-execute::
       :hide-code:

       shape.plot(color='white')
    """

    def __init__(self, *args, **kwargs) -> None:
        kwargs.setdefault("strut_heights", np.sqrt(2.0) / 2.0)
        super().__init__(*args, **kwargs)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        cube_vertices = list(
            product([-self._UNIT_CUBE_SIZE / 2, self._UNIT_CUBE_SIZE / 2], repeat=3)
        )

        face_centers = [
            [sign * self._UNIT_CUBE_SIZE / 2 if i == axis else 0.0 for i in range(3)]
            for axis in range(3)
            for sign in [-1, 1]
        ]

        return np.array(cube_vertices + face_centers)

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        """Generate index pairs representing the struts in the octet-truss using KDTree."""
        tree = KDTree(self.base_vertices)
        pairs = set()
        tolerance = 1e-5

        connection_distance = (self._UNIT_CUBE_SIZE / np.sqrt(2)) + tolerance

        for i, vertex in enumerate(self.base_vertices):
            indices = tree.query_ball_point(vertex, connection_distance)
            for j in indices:
                if i != j:
                    pairs.add(tuple(sorted((i, j))))

        return np.array(list(pairs))
