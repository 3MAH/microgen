"""
Cubic mesh for FE
"""

import warnings
from typing import Optional, Union

import numpy as np
import numpy.typing as npt
import pyvista as pv


class SingleMesh:
    """
    SingleMesh class to manage list of Nodes and Elements inside a Phase
    :param nodes_coords: list of nodes (np.ndarray)
    :param elements : dictionary of elements (key: int, values : np.ndarray). The key is the element type
    :param pvmesh : The pyvista mesh, if it exists already (this could be further generated)
    :param nodes_indices : index of node list (if different from the natural index of nodes array)
    """

    def __init__(
        self,
        nodes_coords: np.ndarray,
        elements: dict[pv.CellType, npt.NDArray[np.int_]],
        pvmesh: Optional[pv.UnstructuredGrid] = None,
        nodes_indices: Optional[npt.NDArray[np.int_]] = None,
    ) -> None:
        self.nodes_coords = nodes_coords
        self.elements = elements  # element dictionary
        self._pvmesh = pvmesh
        self.nodes_indices = (
            nodes_indices  # indices of nodes
        )
        self._surface: Union[None, pv.PolyData] = None

    def _to_cells_and_celltype(self) -> tuple[npt.NDArray[np.int_], npt.NDArray[np.int_]]:
        """
        Returns a numpy array, with the indices of the cells
        The “padding” indicating the number of points per cell is introduced
        exactly as a pyvista UnstructuredGrid.cells property

        :return cells and cells_type : A tuple of flattened arrays that contains all cells and cells_type

        Note that if only tetrahedral elements are handled py SingleMesh, this function
        is able to handle multiple cells for future updates
        """

        cells = np.empty((0,), dtype=int)
        cells_type = np.empty((0,), dtype=int)

        for key in self.elements.keys():
            padding = self.elements[key].shape[1]
            cells_temp = np.empty((self.elements[key].shape[0], padding + 1), dtype=int)
            cells_temp[:, 0] = padding
            cells_temp[:, 1:] = self.elements[key][:, :padding]
            cells_type_temp = np.full((self.elements[key].shape[0]), key)
            cells = np.append(cells, cells_temp)
            cells_type = np.append(cells_type, cells_type_temp)

        return (cells.ravel(), cells_type.ravel())

    def to_pyvista(self) -> pv.UnstructuredGrid:
        """
        Builds a pyvista UnstructuredGrid.
        Node and element data are not copied (by default a shallow copy is operated)
        Note that for now singleMesh works only considering tetrahedral elements

        :return pvmesh : A pyvista.UnstructuredGrid object that contains the (points) nodes,
        cells and celltypes from the SingleMesh object
        """

        cells, celltypes = self._to_cells_and_celltype()

        return pv.UnstructuredGrid(cells, celltypes, self.nodes_coords)

    @staticmethod
    def from_pyvista(pvmesh: Union[pv.UnstructuredGrid, pv.PolyData]):
        """Build a SingleMesh from a pyvista UnstructuredGrid or PolyData mesh (in this last case,
        the PolyData object is cast into an Unstructured Grid).
        Node and element data are not copied (by default a shallow copy is operated)
        Mesh with multiple element types are not handled.
        Note that for now SingleMesh works only considering tetrahedral elements

        :param pvmesh: the mesh as a pyvista UnstructuredGrid object

        :return : A SingleMesh Object from the tetrahedral elements of the pyvista Unstructured Grid
        """
        if isinstance(pvmesh, pv.PolyData):
            pvmesh = pvmesh.cast_to_unstructured_grid()

        try:
            check_if_only_linear_tetrahedral(pvmesh)
            elements = {pv.CellType.TETRA: pvmesh.cells_dict[pv.CellType.TETRA]}
            return SingleMesh(pvmesh.points, elements, pvmesh)
        except ValueError as e:
            print(e)
        return None

    @staticmethod
    def read(filename: str):
        """Build a SingleMesh from a pyvista mesh file. This function uses the pyvista load method.
        The file type is inferred from the file name.

        :param filename : the name of the file that contains the pyvista mesh

        :return : a SingleMesh object
        """

        mesh = SingleMesh.from_pyvista(pv.read(filename))
        return mesh

    @property
    def pvmesh(
        self,
    ) -> pv.UnstructuredGrid:
        """
        Return the pyvista mesh (UnstructuredGrid) of the considered SingleMesh
        """

        if not isinstance(self._pvmesh, pv.UnstructuredGrid):
            self._pvmesh = self.to_pyvista()

        return self._pvmesh

    @property
    def surface(
        self,
    ) -> pv.PolyData:
        """
        Return the surface mesh of the considered mesh
        If it does not exist, the surface mesh is generated using the extract_surface() method from pyvista :
        https://docs.pyvista.org/version/stable/api/core/_autosummary/pyvista.DataSetFilters.extract_surface.html#pyvista.DataSetFilters.extract_surface
        """
        if not isinstance(self._surface, pv.PolyData):
            self._surface = self._extract_surface()

        return self._surface

    def _extract_surface(
        self,
    ) -> pv.PolyData:
        """
        extract the surface mesh of a pv.UnstructuredGrid (stored in self.mesh) using the pyvista extract_surface filter.
        If the mesh as a pyvista Unstructured Grid does not exist, the to_pyvista() method is utilized to generate it

        :return pv.PolyData: surface mesh
        """
        if not (isinstance(self._pvmesh, pv.UnstructuredGrid)):
            self._pvmesh = self.to_pyvista()

        return self._pvmesh.extract_surface()


def check_if_only_linear_tetrahedral(pvmesh: pv.UnstructuredGrid) -> None:
    """
    Check if only linear tetrahedral elements with 10 nodes are present in the pyvista UnstructuredGrid object
    Warnings are prompted if 1D or 2D elements are present since they are not considered
    An error is raised if the mesh contains other than linear_tetrahedral 3D elements

    :param pvmesh: The pyvista UnstructuredGrid object
    """

    # 3:'lin2',
    # 5:'tri3',
    # 9:'quad4',
    # 10:'tet4',
    # 12:'hex8',
    # 13:'wed6',
    # 14:'pyr5',
    # 21:'lin3',
    # 22:'tri6',
    # 23:'quad8',
    # 24:'tet10',
    # 25:'hex20'

    set_elm1d_type = {pv.CellType.LINE, pv.CellType.QUADRATIC_EDGE}
    set_elm2d_type = {pv.CellType.TRIANGLE, pv.CellType.QUAD, pv.CellType.QUADRATIC_TRIANGLE, pv.CellType.QUADRATIC_QUAD}
    set_elm3d_type_other_than_linear_tetra = {pv.CellType.HEXAHEDRON, pv.CellType.WEDGE, pv.CellType.PYRAMID, pv.CellType.QUADRATIC_TETRA, pv.CellType.QUADRATIC_HEXAHEDRON}

    set_cells_in_pvmesh = set(list(pvmesh.cells_dict))

    if set_cells_in_pvmesh.intersection(set_elm1d_type):
        warnings.warn(
            "1D elements are present in the PyVista UnstructuredGrid. They will be ignored."
        )
    if set_cells_in_pvmesh.intersection(set_elm2d_type):
        warnings.warn(
            "2D elements are present in the PyVista UnstructuredGrid. They will be ignored. \nSurface elements shall be extracted automatically from the 3d mesh"
        )

    if set_cells_in_pvmesh.intersection(set_elm3d_type_other_than_linear_tetra):
        raise ValueError(
            "Mesh contains elements other than linear tetrahedra. \nUse triangulate() to ensure that only linear tetrahedra are used."
        )
