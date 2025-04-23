"""Rhombic Dodecahedron.

===============================================================================
Rhombic Dodecahedron (:mod:`microgen.shape.strut_lattice.rhombic_dodecahedron`)
===============================================================================
"""

import math as m
from itertools import product

import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree

from .abstract_lattice import AbstractLattice


class RhombicDodecahedron(AbstractLattice):
    """
    Class to create a unit rhombic dodecahedron lattice of given cell size and density or strut radius

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.RhombicDodecahedron().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, strut_heights=m.sqrt(3.0) / 4.0)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        outer_cube_vertices = list(
            product([-self._UNIT_CUBE_SIZE / 2, self._UNIT_CUBE_SIZE / 2], repeat=3)
        )
        inner_cube_vertices = list(
            product([-self._UNIT_CUBE_SIZE / 4, self._UNIT_CUBE_SIZE / 4], repeat=3)
        )

        face_centers = [
            [sign * self._UNIT_CUBE_SIZE / 2 if i == axis else 0.0 for i in range(3)]
            for axis in range(3)
            for sign in [-1, 1]
        ]

        return np.array(outer_cube_vertices + inner_cube_vertices + face_centers)

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        tree = KDTree(self.base_vertices)
        tolerance = 1e-5
        target_distance = m.sqrt(3.0) / 4.0 + tolerance
        pairs = set(tree.query_pairs(r=target_distance))
        return np.array(list(pairs))
