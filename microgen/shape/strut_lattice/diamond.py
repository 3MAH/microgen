"""
=====================================================
Diamond (:mod:`microgen.shape.strut_lattice.diamond`)
=====================================================
"""

from itertools import product

import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree

from .abstract_lattice import AbstractLattice


class Diamond(AbstractLattice):
    """
    Class to create a unit diamond lattice of given cell size and density or strut radius

    .. jupyter-execute::
       :hide-code:
       :hide-output:

       import microgen

       shape = microgen.Diamond(strut_radius=0.1).generate_vtk()

    .. jupyter-execute::
       :hide-code:

       shape.plot(color='white')
    """

    def __init__(self, *args, **kwargs) -> None:
        self._tetrahedra_centers = self._generate_tetrahedra_centers()
        self._tetrahedra_vertices = self._generate_tetrahedra_vertices()
        kwargs.setdefault("strut_heights", np.sqrt(3.0) / 4.0)
        super().__init__(*args, **kwargs)

    def _generate_tetrahedra_centers(self) -> npt.NDArray[np.float64]:
        candidates = np.array(
            list(
                product([-self._UNIT_CUBE_SIZE / 4, self._UNIT_CUBE_SIZE / 4], repeat=3)
            )
        )
        centers = candidates[np.sum(candidates < 0, axis=1) % 2 == 0]

        return centers

    def _generate_tetrahedra_vertices(self) -> npt.NDArray[np.float64]:
        face_centers = [
            [sign * self._UNIT_CUBE_SIZE / 2 if i == axis else 0.0 for i in range(3)]
            for axis in range(3)
            for sign in [-1, 1]
        ]
        candidates = np.array(
            list(
                product([-self._UNIT_CUBE_SIZE / 2, self._UNIT_CUBE_SIZE / 2], repeat=3)
            )
        )
        outer_cube_corners = candidates[np.sum(candidates < 0, axis=1) % 2 == 0]

        return np.vstack((face_centers, outer_cube_corners))

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        return np.vstack((self._tetrahedra_centers, self._tetrahedra_vertices))

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        tree = KDTree(self.base_vertices)
        pairs = set()
        for i, center in enumerate(self._tetrahedra_centers):
            _, indices = tree.query(center, k=5)
            for j in indices[1:5]:
                pair = tuple(sorted([i, j]))
                pairs.add(pair)

        return np.array(list(pairs))
