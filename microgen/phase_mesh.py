
"""
Cubic mesh for FE
"""

import pyvista as pv
import numpy as np
from .rve import Rve

class PhaseMesh:
    """
    PhaseMesh class to manage list of Nodes and Elements inside a Phase
    :param nodes: list of nodes (np.ndarray)
    :param elements : list of elements (np.ndarray)
    :param elm_type : type of elements (np.ndarray)
    :param nodes_index : index of node list (if different from the natural index of nodes array)       
    """
    def __init__(
        self,
        nodes : Optional[np.ndarray] = None,
        elements : Optional[np.ndarray] = None,
        elm_type : Optional[str] = None,
        nodes_index : Optional[np.ndarray] = None,
    ) -> None:
    
        self.nodes = nodes #node coordinates
        self.elements = elements #element table
        self.elm_type = elm_type
        self.nodes_index = nodes_index

        self.pv_mesh = pv.UnstructuredGrid(self.nodes, self.elements, self.elm_type)