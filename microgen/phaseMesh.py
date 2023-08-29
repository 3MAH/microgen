
"""
Cubic mesh for FE
"""

from typing import Sequence, Optional, Union

import pyvista as pv
import numpy as np
import warnings
from .rve import Rve

class PhaseMesh:
    """
    PhaseMesh class to manage list of Nodes and Elements inside a Phase
    :param nodes: list of nodes (np.ndarray)
    :param elements : dictionnary of elements (key: int, values : np.ndarray). The key is the element type
    :param mesh : The pyvista mesh, if it exists already (this could be further generated)
    :param nodes_index : index of node list (if different from the natural index of nodes array)       
    """
    def __init__(
        self,
        nodes : np.ndarray,
        elements : dict,
        mesh :Optional[pv.UnstructuredGrid] = None,
        name : Optional[str] = None,
        nodes_index : Optional[np.ndarray] = None,
    ) -> None:
    
        self.nodes = nodes #node coordinates
        self.elements = elements #element dictionnary
        self._mesh = mesh
        self.name = name
        self.nodes_index = nodes_index #indices of nodes (i.e, if they come from a bigger mesh)
        self._surface = None

    def _to_cells_and_celltype(self):
        """Return a numpy array, with the indices of the cells
        The “padding” indicating the number of points per cell in introduced
        exactly as a pyvista UnstructuredGrid.cells property"""
    
        cells = np.empty((0,), dtype=int)
        cells_type = np.empty((0,), dtype=int)

        for key in self.elements.keys():
            
            padding = self.elements[key].shape[1]
            cells_temp = np.empty((self.elements[key].shape[0], padding+1), dtype=int)
            cells_temp[:,0] = padding
            cells_temp[:,1:] = self.elements[key][:,:padding]
            cells_type_temp = np.full((self.elements[key].shape[0]), key)
            cells = np.append(cells, cells_temp)
            cells_type = np.append(cells_type, cells_type_temp)            

        return (cells.ravel(), cells_type.ravel())
    
    def to_pyvista(self):
        """Build a pyvista UnstructuredGrid. 
        Node and element data are not copied (to be verified with pyvista).
        Mesh with multiple element type are handled."""

        cells, celltypes = self._to_cells_and_celltype()      

        return pv.UnstructuredGrid(cells, celltypes, self.nodes)

    @staticmethod
    def from_pyvista(pvmesh, name = ""):
        """Build a Mesh from a pyvista UnstructuredGrid or PolyData mesh. 
        Node and element data are not copied.
        Mesh with multiple element type are handled."""                       

        if isinstance(pvmesh, pv.PolyData):
            pvmesh = pvmesh.cast_to_unstructured_grid()
                
        try: 
            _check_if_only_tetra(pvmesh)
            elements = {10: pvmesh.cells_dict(10)}
            return PhaseMesh(pvmesh.points, elements, pvmesh)            
        except ValueError as e:     
            print(e)  

        #extract only the tetrahedral elements
        elements = {10: pvmesh.cells_dict(10)}

    @staticmethod
    def read(
        filename : str,
        name : Optional[str] = None
    ):
        """Build a Mesh from a mesh. This function use the pyvista load method. 
        The file type is inferred from the file name."""

        mesh = PhaseMesh.from_pyvista(pv.read(filename), name=name)
        return mesh

    @property
    def mesh(
        self,
    ) -> None:
        """
        Return the pyvista mesh (UnstructuredGrid) of the considered boxMesh
        """  
        if not(isinstance(self._mesh, pv.UnstructuredGrid)):
            self._mesh = self.to_pyvista()

        return self._mesh

    @mesh.setter
    def mesh(
        self,
        value: pv.UnstructuredGrid,
    ) -> None:
        """
        set a the pyvista mesh

        :param value: Pyvista mesh (UnstructuredGrid)
        """            

        self._mesh = value

    @property
    def surface(
        self,
    ) -> None:
        """
        Return the surface mesh of the considered mesh
        """  
        if not(isinstance(self._surface, pv.PolyData)):       
            self._surface = self._extract_surface()

        return self._surface

    @surface.setter
    def surface(
        self,
        value: pv.PolyData,
    ) -> None:
        """
        set a the surface of the mesh

        :param value: surface of the mesh
        """            

        self._surface = value

    def _extract_surface(
        self,
    ) -> pv.PolyData:
        """
        extract the surface mesh of a pv.UnstructuredGrid using the pyvista extract_surface filter.

        :return pv.PolyData: surface mesh
        """    
        if not(isinstance(self._mesh, pv.UnstructuredGrid)):
            self._mesh = self.to_pyvista()            

        return self._mesh.extract_surface()
    
def _check_if_only_tetra(
    pvmesh : pv.UnstructuredGrid
)-> None:
    
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
    
    set_elm1d_type = {3, 21}
    set_elm2d_type = {5, 9, 22}
    set_elm3d_type_other_than_10 = {10, 12, 13, 14, 23, 24, 25}

    set_cells_in_pvmesh = set(list(pvmesh.cells_dict))

    if set_cells_in_pvmesh.intersection(set_elm1d_type) > 0:
        warnings.warn("1D elements are present in the PyVista UnstructuredGrid. There will be ignored.")
    if set_cells_in_pvmesh.intersection(set_elm2d_type) > 0:
        warnings.warn("2D elements are present in the PyVista UnstructuredGrid. There will be ignored \n", 
            "Surface elements shall be extracted automatically from the 3d mesh")
        
    if set_cells_in_pvmesh.intersection(set_elm3d_type_other_than_10) > 0:
        raise ValueError("Mesh contains elements other than linear tetrahedra \n",
            "You shall use triangulate to ensure that only linear tetrahedra are used ")