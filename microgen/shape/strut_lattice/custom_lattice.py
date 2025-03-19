from .abstract_lattice import AbstractLattice
import numpy as np
import numpy.typing as npt


class CustomLattice(AbstractLattice):
    """
    Class to create a custom lattice with user-defined base vertices and strut vertex pairs
    """

    def __init__(
        self,
        base_vertices: npt.NDArray[np.float64],
        strut_vertex_pairs: npt.NDArray[np.int64],
        *args,
        **kwargs,
    ) -> None:
        super().__init__(
            *args,
            **kwargs,
            base_vertices=base_vertices,
            strut_vertex_pairs=strut_vertex_pairs,
        )

    def _generate_base_vertices(self) -> npt.NDArray[np.float64]:
        return self.base_vertices

    def _generate_strut_vertex_pairs(self) -> npt.NDArray[np.int64]:
        return self.strut_vertex_pairs
