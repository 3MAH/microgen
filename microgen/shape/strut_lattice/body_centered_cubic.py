"""Body Centered Cubic.

=============================================================================
Body Centered Cubic (:mod:`microgen.shape.strut_lattice.body_centered_cubic`)
=============================================================================
"""

from itertools import product

import numpy as np
import numpy.typing as npt

from .abstract_lattice import AbstractLattice


class BodyCenteredCubic(AbstractLattice):
    """
    Class to create a unit body-centered cubic lattice of given cell size and density or strut radius

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.BodyCenteredCubic().generate_vtk()
       shape.plot(color='white')
    """

    def __init__(self, *args, **kwargs) -> None:
        kwargs.setdefault("strut_heights", np.sqrt(3.0) / 2.0)
        super().__init__(*args, **kwargs)

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        unit_cube_vertices = np.array(
            list(
                product([-self._UNIT_CUBE_SIZE / 2, self._UNIT_CUBE_SIZE / 2], repeat=3)
            )
        )
        return np.vstack(([0.0, 0.0, 0.0], unit_cube_vertices))

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        return np.array([[0, i] for i in range(1, 9)])
