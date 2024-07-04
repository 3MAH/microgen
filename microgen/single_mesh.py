"""Cubic mesh for FE."""

from __future__ import annotations

import warnings

import numpy as np
import numpy.typing as npt
import pyvista as pv


class NotOnlyLinearTetrahedraError(Exception):
    """Raised when 3d elements other than linear tetrahedra are found in a mesh."""


class SingleMesh:
    """SingleMesh class to manage list of Nodes and Elements inside a Phase.

    :param nodes_coords: list of nodes (np.ndarray)
    :param elements : dictionary of elements (key: int, values : np.ndarray).
        The key is the element type
    :param nodes_indices : index of node list
        (if different from the natural index of nodes array)
    """

    def __init__(
        self: SingleMesh,
        nodes_coords: npt.NDArray[np.float64],
        elements: dict[pv.CellType, npt.NDArray[np.int_]],
        nodes_indices: npt.NDArray[np.int_] | None = None,
    ) -> None:
        """Initialize the SingleMesh object with nodes and elements."""
        self.nodes_coords = nodes_coords
        self.elements = elements  # element dictionary
        self._pvmesh: pv.UnstructuredGrid | None = None
        self.nodes_indices = nodes_indices  # indices of nodes
        self._surface: pv.PolyData | None = None

    def _to_cells_and_celltype(
        self: SingleMesh,
    ) -> tuple[npt.NDArray[np.int_], npt.NDArray[np.int_]]:
        """Return a numpy array, with the indices of the cells.

        The “padding” indicating the number of points per cell is introduced
        exactly as a pyvista UnstructuredGrid.cells property

        :return cells and cells_type:
            A tuple of flattened arrays that contains all cells and cells_type

        Note that if only tetrahedral elements are handled py SingleMesh, this function
        is able to handle multiple cells for future updates
        """
        cells = np.empty((0,), dtype=int)
        cells_type = np.empty((0,), dtype=int)

        for key in self.elements:
            padding = self.elements[key].shape[1]
            cells_temp = np.empty((self.elements[key].shape[0], padding + 1), dtype=int)
            cells_temp[:, 0] = padding
            cells_temp[:, 1:] = self.elements[key][:, :padding]
            cells_type_temp = np.full((self.elements[key].shape[0]), key)
            cells = np.append(cells, cells_temp)
            cells_type = np.append(cells_type, cells_type_temp)

        return (cells.ravel(), cells_type.ravel())

    def to_pyvista(self: SingleMesh) -> pv.UnstructuredGrid:
        """Build a pyvista UnstructuredGrid.

        Node and element data are not copied (by default a shallow copy is operated)
        Note that for now singleMesh works only considering tetrahedral elements

        :return pvmesh :
            A pyvista.UnstructuredGrid object that contains the (points) nodes,
            cells and celltypes from the SingleMesh object
        """
        cells, celltypes = self._to_cells_and_celltype()

        return pv.UnstructuredGrid(cells, celltypes, self.nodes_coords)

    @staticmethod
    def from_pyvista(pvmesh: pv.UnstructuredGrid | pv.PolyData) -> SingleMesh:
        """Build a SingleMesh from a pyvista UnstructuredGrid or PolyData mesh.

        In the case of a PolyData object, it is cast into an Unstructured Grid).
        Node and element data are not copied (by default a shallow copy is operated)
        Mesh with multiple element types are not handled.
        Note that for now SingleMesh works only considering tetrahedral elements

        :param pvmesh: the mesh as a pyvista UnstructuredGrid object

        :return : A SingleMesh Object from the tetrahedral elements of
            the pyvista Unstructured Grid
        """
        if isinstance(pvmesh, pv.PolyData):
            pvmesh = pvmesh.cast_to_unstructured_grid()

        # extract only the tetrahedral elements
        # raises an Exception if 3d elements other than
        # linear tetrahedra are found in the mesh
        check_if_only_linear_tetrahedral(pvmesh)
        elements = {pv.CellType.TETRA: pvmesh.cells_dict[pv.CellType.TETRA]}
        return SingleMesh(pvmesh.points, elements)

    @staticmethod
    def read(filename: str) -> SingleMesh:
        """Build a SingleMesh from a pyvista mesh file.

        This function uses the pyvista load method.
        The file type is inferred from the file name.

        :param filename : the name of the file that contains the pyvista mesh

        :return : a SingleMesh object
        """
        return SingleMesh.from_pyvista(pv.read(filename))

    @property
    def pvmesh(self: SingleMesh) -> pv.UnstructuredGrid:
        """Return the pyvista mesh (UnstructuredGrid) of the considered SingleMesh."""
        if not isinstance(self._pvmesh, pv.UnstructuredGrid):
            self._pvmesh = self.to_pyvista()

        return self._pvmesh

    @property
    def surface(self: SingleMesh) -> pv.PolyData:
        """Return the surface mesh of the considered mesh.

        If it does not exist, the surface mesh is generated using the extract_surface()
        method from pyvista : https://docs.pyvista.org/version/stable/api/core/_autosummary/pyvista.DataSetFilters.extract_surface.html#pyvista.DataSetFilters.extract_surface
        """
        if not isinstance(self._surface, pv.PolyData):
            self._surface = self._extract_surface()

        return self._surface

    def _extract_surface(
        self: SingleMesh,
    ) -> pv.PolyData:
        """Extract the surface mesh of a pv.UnstructuredGrid.

        Stored in self.mesh using the pyvista extract_surface filter
        If the mesh as a pyvista Unstructured Grid does not exist, the to_pyvista()
        method is utilized to generate it

        :return pv.PolyData: surface mesh
        """
        if not isinstance(self._pvmesh, pv.UnstructuredGrid):
            self._pvmesh = self.to_pyvista()

        return self._pvmesh.extract_surface()


def check_if_only_linear_tetrahedral(pvmesh: pv.UnstructuredGrid) -> None:
    """Check if only linear tetrahedral elements in the mesh.

    Check if only tetra elements with 10 nodes are present in the
    pyvista UnstructuredGrid object.
    Warnings are prompted if 1D or 2D elements are present since they are not considered
    An error is raised if the mesh contains other than linear_tetrahedral 3D elements

    :param pvmesh: The pyvista UnstructuredGrid object
    """
    set_elm1d_type = {pv.CellType.LINE, pv.CellType.QUADRATIC_EDGE}
    set_elm2d_type = {
        pv.CellType.TRIANGLE,
        pv.CellType.QUAD,
        pv.CellType.QUADRATIC_TRIANGLE,
        pv.CellType.QUADRATIC_QUAD,
    }
    set_elm3d_type_other_than_linear_tetra = {
        pv.CellType.HEXAHEDRON,
        pv.CellType.WEDGE,
        pv.CellType.PYRAMID,
        pv.CellType.QUADRATIC_TETRA,
        pv.CellType.QUADRATIC_HEXAHEDRON,
    }

    set_cells_in_pvmesh = set(pvmesh.cells_dict)

    if set_cells_in_pvmesh.intersection(set_elm1d_type):
        warnings.warn(
            (
                "1D elements are present in the PyVista UnstructuredGrid. "
                "They will be ignored."
            ),
            stacklevel=2,
        )
    if set_cells_in_pvmesh.intersection(set_elm2d_type):
        warnings.warn(
            (
                "2D elements are present in the PyVista UnstructuredGrid. "
                "They will be ignored. "
                "\nSurface elements shall be extracted automatically from the 3d mesh."
            ),
            stacklevel=2,
        )

    if set_cells_in_pvmesh.intersection(set_elm3d_type_other_than_linear_tetra):
        err_msg = (
            "Mesh contains elements other than linear tetrahedra."
            "Use triangulate() to ensure that only linear tetrahedra are used."
        )
        raise NotOnlyLinearTetrahedraError(err_msg)
