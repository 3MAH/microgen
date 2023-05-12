
"""
Cubic mesh for FE
"""

from typing import Sequence, Optional, Union

import pyvista as pv
import numpy as np
from .rve import Rve

class PhaseMesh:
    """
    PhaseMesh class to manage list of Nodes and Elements inside a Phase
    :param nodes: list of nodes (np.ndarray)
    :param elements : dictionnary of elements (key: int, values : np.ndarray). The key is the element type
    :param nodes_index : index of node list (if different from the natural index of nodes array)       
    """
    def __init__(
        self,
        nodes : np.ndarray,
        elements : dict,
        name : Optional[np.ndarray] = None,
        nodes_index : Optional[np.ndarray] = None,
    ) -> None:
    
        self.nodes = nodes #node coordinates
        self.elements = elements #element dictionnary
        self.name = name
        self.nodes_index = nodes_index #indices of nodes (i.e, if they come from a bigger mesh)

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
            print(cells_temp)
            print(cells_type_temp)    
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
        
        elements = pvmesh.cells_dict
        return PhaseMesh(pvmesh.points, elements)
    
    @staticmethod
    def read(filename, name = ""):
        """Build a Mesh from a mesh. This function use the pyvista load method. 
        The file type is inferred from the file name."""

        mesh = PhaseMesh.from_pyvista(pv.read(filename), name=name)
        return mesh
