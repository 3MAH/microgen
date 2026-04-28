"""
=================================================================
Cuboctahedron (:mod:`microgen.shape.strut_lattice.cuboctahedron`)
=================================================================
"""

from itertools import product

import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree

from .abstract_lattice import BALL_POINT_RADIUS_TOLERANCE, AbstractLattice


class Cuboctahedron(AbstractLattice):
    """
    Class to create a unit cuboctahedron lattice of given cell size and density or strut radius

    .. jupyter-execute::
       :hide-code:
       :hide-output:

       import microgen

       shape = microgen.Cuboctahedron(strut_radius=0.1).generate_vtk()

    .. jupyter-execute::
       :hide-code:

       shape.plot(color='white')
    """

    def __init__(self, *args, **kwargs) -> None:
        kwargs.setdefault("strut_heights", np.sqrt(2.0) / 2.0)
        super().__init__(*args, **kwargs)

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
        for i, indices in enumerate(
            tree.query_ball_point(
                self.base_vertices,
                r=self._UNIT_CUBE_SIZE / np.sqrt(2.0) + BALL_POINT_RADIUS_TOLERANCE,
            )
        ):
            for j in indices:
                if i != j:
                    pairs.add(tuple(sorted((i, j))))

        return np.array(list(pairs))
