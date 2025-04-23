"""Cuboctahedron.

=================================================================
Cuboctahedron (:mod:`microgen.shape.strut_lattice.cuboctahedron`)
=================================================================
"""

import math as m
from itertools import product

import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree

from .abstract_lattice import AbstractLattice


class Cuboctahedron(AbstractLattice):
    """
    Class to create a unit cuboctahedron lattice of given cell size and density or strut radius

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Cuboctahedron().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, strut_heights=m.sqrt(2.0) / 2.0)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        cube_vertices = np.array(
            list(
                product([-self._UNIT_CUBE_SIZE / 2, self._UNIT_CUBE_SIZE / 2], repeat=3)
            )
        )

        edges = [
            (i, j)
            for i in range(len(cube_vertices))
            for j in range(i + 1, len(cube_vertices))
            if np.sum(np.abs(cube_vertices[i] - cube_vertices[j]))
            == self._UNIT_CUBE_SIZE
        ]

        return np.array([(cube_vertices[i] + cube_vertices[j]) / 2.0 for i, j in edges])

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        tree = KDTree(self.base_vertices)
        pairs = set()
        tolerance = 1e-5
        for i, indices in enumerate(
            tree.query_ball_point(
                self.base_vertices, r=self._UNIT_CUBE_SIZE / m.sqrt(2.0) + tolerance
            )
        ):
            for j in indices:
                if i != j:
                    pairs.add(tuple(sorted((i, j))))

        return np.array(list(pairs))
