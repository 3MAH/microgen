
"""
Cubic mesh for FE
"""

from typing import Sequence, Optional, Union

import pyvista as pv
import numpy as np
from .rve import Rve
from .phaseMesh import PhaseMesh

class CubicMesh(PhaseMesh):
    """
    CubicMesh class to manage list of Nodes and Elements inside an Rve
    :param node_list: list of nodes
    :param node_list_name : name of node lists
    """
    def __init__(
        self,
        nodes : np.ndarray,
        elements : np.ndarray,
        elm_type : np.ndarray,
        nodes_index : Optional[np.ndarray] = None,
    ) -> None:
        super().__init__(nodes=nodes, elements=elements, elm_type=elm_type, nodes_index=nodes_index)

        self.center = None

        self.corner_list_XmYmZm = None
        self.corner_list_XmYpZm = None
        self.corner_list_XpYmZm = None
        self.corner_list_XpYpZm = None
        self.corner_list_XmYmZp = None
        self.corner_list_XmYpZp = None
        self.corner_list_XpYmZp = None
        self.corner_list_XpYpZp = None

        self.corners = []

        self.edge_list_XmYm = None
        self.edge_list_XpYm = None
        self.edge_list_XpYp = None
        self.edge_list_XmYp = None
        self.edge_list_XmZm = None
        self.edge_list_XpZm = None
        self.edge_list_XpZp = None
        self.edge_list_XmZp = None
        self.edge_list_YmZm = None
        self.edge_list_YpZm = None
        self.edge_list_YpZp = None
        self.edge_list_YmZp = None

        self.edges = []

        self.face_list_Xm = None
        self.face_list_Ym = None
        self.face_list_Zm = None
        self.face_list_Xp = None
        self.face_list_Yp = None
        self.face_list_Zp = None

        self.faces = []


    def construct(
        self,
        rve: Rve,
        tol: float = 1.e-8,
    ) -> None:
            
        crd = self.nodes
        pv_mesh = self.pv_mesh
        closest_point_to_rve_center = pv_mesh.find_closest_point(rve.center)
        self.center = crd[closest_point_to_rve_center]

        self.face_list_Xm = np.where( np.abs(crd[:,0] - rve.x_min) < tol )[0]
        self.face_list_Xp = np.where( np.abs(crd[:,0] - rve.x_max) < tol )[0]
        
        self.face_list_Ym = np.where( np.abs(crd[:,1] - rve.y_min) < tol )[0]
        self.face_list_Yp = np.where( np.abs(crd[:,1] - rve.y_max) < tol )[0]

        self.face_list_Zm = np.where( np.abs(crd[:,2] - rve.z_min) < tol )[0]
        self.face_list_Zp = np.where( np.abs(crd[:,2] - rve.z_max) < tol )[0]

        self.edge_list_XmYm = np.intersect1d(self.face_list_Xm , self.face_list_Ym, assume_unique=True)
        self.edge_list_XpYm = np.intersect1d(self.face_list_Xp , self.face_list_Ym, assume_unique=True)
        self.edge_list_XpYp = np.intersect1d(self.face_list_Xp , self.face_list_Yp, assume_unique=True)
        self.edge_list_XmYp = np.intersect1d(self.face_list_Xm , self.face_list_Yp, assume_unique=True)
        
        self.edge_list_XmZm = np.intersect1d(self.face_list_Xm , self.face_list_Zm, assume_unique=True)
        self.edge_list_XpZm = np.intersect1d(self.face_list_Xp , self.face_list_Zm, assume_unique=True)
        self.edge_list_XpZp = np.intersect1d(self.face_list_Xp , self.face_list_Zp, assume_unique=True)
        self.edge_list_XmZp = np.intersect1d(self.face_list_Xm , self.face_list_Zp, assume_unique=True)
        
        self.edge_list_YmZm = np.intersect1d(self.face_list_Ym , self.face_list_Zm, assume_unique=True)
        self.edge_list_YpZm = np.intersect1d(self.face_list_Yp , self.face_list_Zm, assume_unique=True)
        self.edge_list_YpZp = np.intersect1d(self.face_list_Yp , self.face_list_Zp, assume_unique=True)
        self.edge_list_YmZp = np.intersect1d(self.face_list_Ym , self.face_list_Zp, assume_unique=True)

        #extract corners from the intersection of edges
        self.corner_list_XmYmZm = np.intersect1d(self.edge_list_XmYm , self.edge_list_YmZm, assume_unique=True)
        self.corner_list_XmYpZm = np.intersect1d(self.edge_list_XmYp , self.edge_list_YpZm, assume_unique=True)
        self.corner_list_XpYmZm = np.intersect1d(self.edge_list_XpYm , self.edge_list_YmZm, assume_unique=True)
        self.corner_list_XpYpZm = np.intersect1d(self.edge_list_XpYp , self.edge_list_YpZm, assume_unique=True)
        self.corner_list_XmYmZp = np.intersect1d(self.edge_list_XmYm , self.edge_list_YmZp, assume_unique=True)
        self.corner_list_XmYpZp = np.intersect1d(self.edge_list_XmYp , self.edge_list_YpZp, assume_unique=True)
        self.corner_list_XpYmZp = np.intersect1d(self.edge_list_XpYm , self.edge_list_YmZp, assume_unique=True)
        self.corner_list_XpYpZp = np.intersect1d(self.edge_list_XpYp , self.edge_list_YpZp, assume_unique=True)

        # Remove nodes that beloing to several sets
        all_corners = np.hstack((self.corner_list_XmYmZm, self.corner_list_XmYpZm, self.corner_list_XpYmZm, self.corner_list_XpYpZm,
                              self.corner_list_XmYmZp, self.corner_list_XmYpZp, self.corner_list_XpYmZp, self.corner_list_XpYpZp))

        self.edge_list_XmYm = np.setdiff1d(self.edge_list_XmYm, all_corners, assume_unique=True)
        self.edge_list_XpYm = np.setdiff1d(self.edge_list_XpYm, all_corners, assume_unique=True)
        self.edge_list_XpYp = np.setdiff1d(self.edge_list_XpYp, all_corners, assume_unique=True)
        self.edge_list_XmYp = np.setdiff1d(self.edge_list_XmYp, all_corners, assume_unique=True)

        self.edge_list_XmZm = np.setdiff1d(self.edge_list_XmZm, all_corners, assume_unique=True)
        self.edge_list_XpZm = np.setdiff1d(self.edge_list_XpZm, all_corners, assume_unique=True)
        self.edge_list_XpZp = np.setdiff1d(self.edge_list_XpZp, all_corners, assume_unique=True)
        self.edge_list_XmZp = np.setdiff1d(self.edge_list_XmZp, all_corners, assume_unique=True)

        self.edge_list_YmZm = np.setdiff1d(self.edge_list_YmZm, all_corners, assume_unique=True)
        self.edge_list_YpZm = np.setdiff1d(self.edge_list_YpZm, all_corners, assume_unique=True)
        self.edge_list_YpZp = np.setdiff1d(self.edge_list_YpZp, all_corners, assume_unique=True)
        self.edge_list_YmZp = np.setdiff1d(self.edge_list_YmZp, all_corners, assume_unique=True)

        all_edges_corners = np.hstack((self.edge_list_XmYm, self.edge_list_XpYm, self.edge_list_XpYp, self.edge_list_XmYp, self.edge_list_XmZm, self.edge_list_XpZm,
                                self.edge_list_XpZp, self.edge_list_XmZp, self.edge_list_YmZm, self.edge_list_YpZm, self.edge_list_YpZp, self.edge_list_YmZp,
                                all_corners))

        self.face_list_Xm = np.setdiff1d(self.face_list_Xm, all_edges_corners, assume_unique=True)
        self.face_list_Xp = np.setdiff1d(self.face_list_Xp, all_edges_corners, assume_unique=True)
        self.face_list_Ym = np.setdiff1d(self.face_list_Ym, all_edges_corners, assume_unique=True)
        self.face_list_Yp = np.setdiff1d(self.face_list_Yp, all_edges_corners, assume_unique=True)
        self.face_list_Zm = np.setdiff1d(self.face_list_Zm, all_edges_corners, assume_unique=True)
        self.face_list_Zp = np.setdiff1d(self.face_list_Zp, all_edges_corners, assume_unique=True)

        self.edge_list_XmYm = self.edge_list_XmYm[np.argsort(crd[self.edge_list_XmYm,2])]
        self.edge_list_XpYm = self.edge_list_XpYm[np.argsort(crd[self.edge_list_XpYm,2])]
        self.edge_list_XpYp = self.edge_list_XpYp[np.argsort(crd[self.edge_list_XpYp,2])]
        self.edge_list_XmYp = self.edge_list_XmYp[np.argsort(crd[self.edge_list_XmYp,2])]
        
        self.edge_list_XmZm = self.edge_list_XmZm[np.argsort(crd[self.edge_list_XmZm,1])]
        self.edge_list_XpZm = self.edge_list_XpZm[np.argsort(crd[self.edge_list_XpZm,1])]
        self.edge_list_XpZp = self.edge_list_XpZp[np.argsort(crd[self.edge_list_XpZp,1])]
        self.edge_list_XmZp = self.edge_list_XmZp[np.argsort(crd[self.edge_list_XmZp,1])]
        
        self.edge_list_YmZm = self.edge_list_YmZm[np.argsort(crd[self.edge_list_YmZm,0])]
        self.edge_list_YpZm = self.edge_list_YpZm[np.argsort(crd[self.edge_list_YpZm,0])]
        self.edge_list_YpZp = self.edge_list_YpZp[np.argsort(crd[self.edge_list_YpZp,0])]
        self.edge_list_YmZp = self.edge_list_YmZp[np.argsort(crd[self.edge_list_YmZp,0])]  

        decimal_round = int(-np.log10(tol)-1)
        self.face_list_Xm = self.face_list_Xm[np.lexsort((crd[self.face_list_Xm, 1], crd[self.face_list_Xm, 2].round(decimal_round)))]
        self.face_list_Xp = self.face_list_Xp[np.lexsort((crd[self.face_list_Xp, 1], crd[self.face_list_Xp, 2].round(decimal_round)))]
        self.face_list_Ym = self.face_list_Ym[np.lexsort((crd[self.face_list_Ym, 0], crd[self.face_list_Ym, 2].round(decimal_round)))]
        self.face_list_Yp = self.face_list_Yp[np.lexsort((crd[self.face_list_Yp, 0], crd[self.face_list_Yp, 2].round(decimal_round)))]
        self.face_list_Zm = self.face_list_Zm[np.lexsort((crd[self.face_list_Zm, 0], crd[self.face_list_Zm, 1].round(decimal_round)))]
        self.face_list_Zp = self.face_list_Zp[np.lexsort((crd[self.face_list_Zp ,0], crd[self.face_list_Zp ,1].round(decimal_round)))]

        print(self.edge_list_XmYm)
        print("ouille\n")
    
