
"""
Cubic mesh for FE
"""

import pyvista as pv
import numpy as np
from .rve import Rve
from .phase_mesh import PhaseMesh

class CubicMesh(PhaseMesh):
    """
    CubicMesh class to manage list of Nodes and Elements inside an Rve
    :param node_list: list of nodes
    :param node_list_name : name of node lists
    """
    def __init__(
        self,
        nodes : Optional[np.ndarray] = None,
        elements : Optional[np.ndarray] = None,
        elm_type : Optional[str] = None,
        nodes_index : Optional[np.ndarray] = None,
    ) -> None:
        super.__init__(nodes=nodes, center=elements, orientation=elm_type, nodes_index=nodes_index)

        self.center = None

        self.corner_listXmYmZm = None
        self.corner_listXmYpZm = None
        self.corner_listXpYmZm = None
        self.corner_listXpYpZm = None
        self.corner_listXmYmZp = None
        self.corner_listXmYpZp = None
        self.corner_listXpYmZp = None
        self.corner_listXpYpZp = None

        self.corners = []

        self.edge_listXmYm = None
        self.edge_listXpYm = None
        self.edge_listXpYp = None
        self.edge_listXmYp = None
        self.edge_listXmZm = None
        self.edge_listXpZm = None
        self.edge_listXpZp = None
        self.edge_listXmZp = None
        self.edge_listYmZm = None
        self.edge_listYpZm = None
        self.edge_listYpZp = None
        self.edge_listYmZp = None

        self.edges = []

        self.face_listXm = None
        self.face_listYm = None
        self.face_listZm = None
        self.face_listXp = None
        self.face_listYp = None
        self.face_listZp = None

        self.faces = []


    def construct(
        self,
        rve: Rve,
        tol: float = 1.e8,
    ) -> None:
            
        crd = self.nodes
        pv_mesh = self.pv_mesh
        closest_point_to_rve_center = pv_mesh.find_closest_point(rve.center)
        self.center = np.array(closest_point_to_rve_center[:3])    

        self.face_list_Xm = np.where( np.abs(crd[:,0] - xmin) < tol )[0]
        self.face_list_Xp = np.where( np.abs(crd[:,0] - xmax) < tol )[0]
        
        self.face_list_Ym = np.where( np.abs(crd[:,1] - ymin) < tol )[0]
        self.face_list_Yp = np.where( np.abs(crd[:,1] - ymax) < tol )[0]

        self.face_list_Zm = np.where( np.abs(crd[:,2] - zmin) < tol )[0]
        self.face_list_Zp = np.where( np.abs(crd[:,2] - zmax) < tol )[0]

        self.Edge_list_XmYm = np.intersect1d(self.face_list_Xm , self.face_list_Ym, assume_unique=True)
        self.Edge_list_XpYm = np.intersect1d(self.face_list_Xp , self.face_list_Ym, assume_unique=True)
        self.Edge_list_XpYp = np.intersect1d(self.face_list_Xp , self.face_list_Yp, assume_unique=True)
        self.Edge_list_XmYp = np.intersect1d(self.face_list_Xm , self.face_list_Yp, assume_unique=True)
        
        self.Edge_list_XmZm = np.intersect1d(self.face_list_Xm , self.face_list_Zm, assume_unique=True)
        self.Edge_list_XpZm = np.intersect1d(self.face_list_Xp , self.face_list_Zm, assume_unique=True)
        self.Edge_list_XpZp = np.intersect1d(self.face_list_Xp , self.face_list_Zp, assume_unique=True)
        self.Edge_list_XmZp = np.intersect1d(self.face_list_Xm , self.face_list_Zp, assume_unique=True)
        
        self.Edge_list_YmZm = np.intersect1d(self.face_list_Ym , self.face_list_Zm, assume_unique=True)
        self.Edge_list_YpZm = np.intersect1d(self.face_list_Yp , self.face_list_Zm, assume_unique=True)
        self.Edge_list_YpZp = np.intersect1d(self.face_list_Yp , self.face_list_Zp, assume_unique=True)
        self.Edge_list_YmZp = np.intersect1d(self.face_list_Ym , self.face_list_Zp, assume_unique=True)

        #extract corners from the intersection of edges
        self.Corner_listXmYmZm = np.intersect1d(self.Edge_list_XmYm , self.Edge_list_YmZm, assume_unique=True)
        self.Corner_listXmYpZm = np.intersect1d(self.Edge_list_XmYp , self.Edge_list_YpZm, assume_unique=True)
        self.Corner_listXpYmZm = np.intersect1d(self.Edge_list_XpYm , self.Edge_list_YmZm, assume_unique=True)
        self.Corner_listXpYpZm = np.intersect1d(self.Edge_list_XpYp , self.Edge_list_YpZm, assume_unique=True)
        self.Corner_listXmYmZp = np.intersect1d(self.Edge_list_XmYm , self.Edge_list_YmZp, assume_unique=True)
        self.Corner_listXmYpZp = np.intersect1d(self.Edge_list_XmYp , self.Edge_list_YpZp, assume_unique=True)
        self.Corner_listXpYmZp = np.intersect1d(self.Edge_list_XpYm , self.Edge_list_YmZp, assume_unique=True)
        self.Corner_listXpYpZp = np.intersect1d(self.Edge_list_XpYp , self.Edge_list_YpZp, assume_unique=True)

        # Remove nodes that beloing to several sets
        all_corners = np.hstack((self.Corner_listXmYmZm, self.Corner_listXmYpZm, self.Corner_listXpYmZm, self.Corner_listXpYpZm,
                              self.Corner_listXmYmZp, self.Corner_listXmYpZp, self.Corner_listXpYmZp, self.Corner_listXpYpZp))

        self.Edge_list_XmYm = np.setdiff1d(self.Edge_list_XmYm, all_corners, assume_unique=True)
        self.Edge_list_XpYm = np.setdiff1d(self.Edge_list_XpYm, all_corners, assume_unique=True)
        self.Edge_list_XpYp = np.setdiff1d(self.Edge_list_XpYp, all_corners, assume_unique=True)
        self.Edge_list_XmYp = np.setdiff1d(self.Edge_list_XmYp, all_corners, assume_unique=True)

        self.Edge_list_XmZm = np.setdiff1d(self.Edge_list_XmZm, all_corners, assume_unique=True)
        self.Edge_list_XpZm = np.setdiff1d(self.Edge_list_XpZm, all_corners, assume_unique=True)
        self.Edge_list_XpZp = np.setdiff1d(self.Edge_list_XpZp, all_corners, assume_unique=True)
        self.Edge_list_XmZp = np.setdiff1d(self.Edge_list_XmZp, all_corners, assume_unique=True)

        self.Edge_list_YmZm = np.setdiff1d(self.Edge_list_YmZm, all_corners, assume_unique=True)
        self.Edge_list_YpZm = np.setdiff1d(self.Edge_list_YpZm, all_corners, assume_unique=True)
        self.Edge_list_YpZp = np.setdiff1d(self.Edge_list_YpZp, all_corners, assume_unique=True)
        self.Edge_list_YmZp = np.setdiff1d(self.Edge_list_YmZp, all_corners, assume_unique=True)
        
        all_edges_corners = np.hstack((self.Edge_list_XmYm, self.Edge_list_XpYm, self.Edge_list_XpYp, self.Edge_list_XmYp, self.Edge_list_XmZm, self.Edge_list_XpZm,
                                self.Edge_list_XpZp, self.Edge_list_XmZp, self.Edge_list_YmZm, self.Edge_list_YpZm, self.Edge_list_YpZp, self.Edge_list_YmZp,
                                all_corners))

        self.face_list_Xm = np.setdiff1d(self.face_list_Xm, all_edges_corners, assume_unique=True)
        self.face_list_Xp = np.setdiff1d(self.face_list_Xp, all_edges_corners, assume_unique=True)
        self.face_list_Ym = np.setdiff1d(self.face_list_Ym, all_edges_corners, assume_unique=True)
        self.face_list_Yp = np.setdiff1d(self.face_list_Yp, all_edges_corners, assume_unique=True)
        self.face_list_Zm = np.setdiff1d(self.face_list_Zm, all_edges_corners, assume_unique=True)
        self.face_list_Zp = np.setdiff1d(self.face_list_Zp, all_edges_corners, assume_unique=True)

        self.Edge_list_XmYm = self.Edge_list_XmYm[np.argsort(crd[bottom_back,2])]
        self.Edge_list_XpYm = self.Edge_list_XpYm[np.argsort(crd[bottom_back,2])]
        self.Edge_list_XpYp = self.Edge_list_XpYp[np.argsort(crd[bottom_back,2])]
        self.Edge_list_XmYp = self.Edge_list_XmYp[np.argsort(crd[bottom_back,2])]
        
        self.Edge_list_XmZm = self.Edge_list_XmZm[np.argsort(crd[bottom_back,1])]
        self.Edge_list_XpZm = self.Edge_list_XpZm[np.argsort(crd[bottom_back,1])]
        self.Edge_list_XpZp = self.Edge_list_XpZp[np.argsort(crd[bottom_back,1])]
        self.Edge_list_XmZp = self.Edge_list_XmZp[np.argsort(crd[bottom_back,1])]
        
        self.Edge_list_YmZm = self.Edge_list_YmZm[np.argsort(crd[bottom_back,0])]
        self.Edge_list_YpZm = self.Edge_list_YpZm[np.argsort(crd[bottom_back,0])]
        self.Edge_list_YpZp = self.Edge_list_YpZp[np.argsort(crd[bottom_back,0])]
        self.Edge_list_YmZp = self.Edge_list_YmZp[np.argsort(crd[bottom_back,0])]

        decimal_round = int(-np.log10(tol)-1)
        self.face_list_Xm = self.face_list_Xm[np.lexsort((crd[left  ,1], crd[left  ,2].round(decimal_round)))]
        self.face_list_Xp = self.face_list_Xp[np.lexsort((crd[right ,1], crd[right ,2].round(decimal_round)))]
        self.face_list_Ym = self.face_list_Ym[np.lexsort((crd[bottom,0], crd[bottom,2].round(decimal_round)))]
        self.face_list_Yp = self.face_list_Yp[np.lexsort((crd[top   ,0], crd[top   ,2].round(decimal_round)))]
        self.face_list_Zm = self.face_list_Zm[np.lexsort((crd[back  ,0], crd[back  ,1].round(decimal_round)))]
        self.face_list_Zp = self.face_list_Zp[np.lexsort((crd[front ,0], crd[front ,1].round(decimal_round)))]
    
