
"""
Cubic mesh for FE
"""

from typing import Sequence, Optional, Union

import pyvista as pv
import numpy as np
from scipy.spatial import KDTree
from .rve import Rve
from .phaseMesh import PhaseMesh

class BoxMesh(PhaseMesh):
    """
    CubicMesh class to manage list of Nodes and Elements inside an Rve
    :param nodes: list of nodes (np.ndarray)
    :param elements : list of elements (np.ndarray)
    :param elm_type : type of elements (np.ndarray)
    :param nodes_index : index of node list (if different from the natural index of nodes array)   
    """
    def __init__(
        self,
        nodes : np.ndarray,
        elements : dict,
        name : Optional[np.ndarray] = None,
        nodes_index : Optional[np.ndarray] = None,
    ) -> None:
        super().__init__(nodes=nodes, elements=elements, name=name, nodes_index=nodes_index)

        self._rve = None

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
        rve: Optional[Rve] = None,
        tol: Optional[float] = 1.e-8,
    ) -> None:
            
        if(isinstance(rve, Rve) == False) : 
            rve = self.rve

        crd = self.nodes
        closest_point_to_rve_center = np.linalg.norm(crd-rve.center, axis=1).argmin()
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

        self.corners = np.vstack((self.corner_list_XmYmZm, self.corner_list_XmYpZm, self.corner_list_XpYmZm, self.corner_list_XpYpZm,
                              self.corner_list_XmYmZp, self.corner_list_XmYpZp, self.corner_list_XpYmZp, self.corner_list_XpYpZp))  

        self.edges = np.vstack((self.edge_list_XmYm, self.edge_list_XpYm, self.edge_list_XpYp, self.edge_list_XmYp, self.edge_list_XmZm, self.edge_list_XpZm,
                                self.edge_list_XpZp, self.edge_list_XmZp, self.edge_list_YmZm, self.edge_list_YpZm, self.edge_list_YpZp, self.edge_list_YmZp))

        self.faces = np.vstack((self.face_list_Xm, self.face_list_Xp, self.face_list_Ym, self.face_list_Yp, self.face_list_Zm, self.face_list_Zp))

    @property
    def rve(
        self,
    ) -> Rve:
        """Return a representative volume element (Rve)
        of the considered mesh"""

        if(isinstance(self._rve, Rve)):
            return self._rve
        else:
            return self._build_rve()

    @rve.setter
    def rve(
        self,
        value: Rve,
    ) -> None:
        self._rve = value

    def _build_rve(
        self,
    ) -> Rve:
        pvmesh = self.to_pyvista()
        xmin, xmax, ymin, ymax, zmin, zmax = pvmesh.bounds
        return Rve.from_min_max(xmin, xmax, ymin, ymax, zmin, zmax)

    def closest_points_on_faces(
        self,
        k_neighbours: float = 3,
        rve: Optional[Rve] = None,
        tol: Optional[float] = 1.e-8,
    ) -> None:
            
        if(isinstance(rve, Rve) == False) : 
            rve = self.rve       

        crd = self.nodes

        all_face_Xp = np.hstack((self.face_list_Xp,self.edge_list_XpYm,self.edge_list_XpYp,self.edge_list_XpZm,self.edge_list_XpZp))
        all_face_Yp = np.hstack((self.face_list_Yp,self.edge_list_YpZm,self.edge_list_YpZp,self.edge_list_XmYp,self.edge_list_XpYp))
        all_face_Zp = np.hstack((self.face_list_Zp,self.edge_list_YmZp,self.edge_list_YpZp,self.edge_list_XmZp,self.edge_list_XpZp))

        kdtree_Xp = KDTree(crd[all_face_Xp])
        kdtree_Yp = KDTree(crd[all_face_Yp])
        kdtree_Zp = KDTree(crd[all_face_Zp])

        offset = np.array([rve.dim_x, 0.0, 0.0])
        crd_XmXp = crd[self.face_list_Xm] + offset
        offset = np.array([0.0, rve.dim_y, 0.0])        
        crd_YmYp = crd[self.face_list_Ym] + offset   
        offset = np.array([0.0, 0.0, rve.dim_z])  
        crd_ZmZp = crd[self.face_list_Zm] + offset              

        dist_XmXp, index_XmXp = kdtree_Xp.query(crd_XmXp, k_neighbours)
        dist_YmYp, index_YmYp = kdtree_Yp.query(crd_YmYp, k_neighbours)
        dist_ZmZp, index_ZmZp = kdtree_Zp.query(crd_ZmZp, k_neighbours)        

        return (all_face_Xp[index_XmXp], dist_XmXp, all_face_Yp[index_YmYp], dist_YmYp, all_face_Zp[index_ZmZp], dist_ZmZp)

    def closest_points_on_edges(
        self,
        rve: Optional[Rve] = None,
        tol: Optional[float] = 1.e-8,
    ) -> None:
            
        if(isinstance(rve, Rve) == False) : 
            rve = self.rve       

        crd = self.nodes

        all_edge_XpYm = np.hstack((self.edge_list_XpYm,self.corner_list_XpYmZm,self.corner_list_XpYmZp))
        all_edge_XpYp = np.hstack((self.edge_list_XpYp,self.corner_list_XpYpZm,self.corner_list_XpYpZp))
        all_edge_XmYp = np.hstack((self.edge_list_XmYp,self.corner_list_XmYpZm,self.corner_list_XmYpZp))

        all_edge_XpZm = np.hstack((self.edge_list_XpZm,self.corner_list_XpYmZm,self.corner_list_XpYpZm))
        all_edge_XpZp = np.hstack((self.edge_list_XpZp,self.corner_list_XpYmZp,self.corner_list_XpYpZp))                
        all_edge_XmZp = np.hstack((self.edge_list_XmZp,self.corner_list_XmYmZp,self.corner_list_XmYpZp))                

        all_edge_YpZm = np.hstack((self.edge_list_YmZm,self.corner_list_XmYpZm,self.corner_list_XpYpZm))
        all_edge_YpZp = np.hstack((self.edge_list_YpZp,self.corner_list_XmYpZp,self.corner_list_XpYpZp))        
        all_edge_YmZp = np.hstack((self.edge_list_YmZp,self.corner_list_XmYmZp,self.corner_list_XpYmZp))        

        kdtree_XpYm = KDTree(crd[all_edge_XpYm])
        kdtree_XpYp = KDTree(crd[all_edge_XpYp])
        kdtree_XmYp = KDTree(crd[all_edge_XmYp])

        kdtree_XpZm = KDTree(crd[all_edge_XpZm])
        kdtree_XpZp = KDTree(crd[all_edge_XpZp])
        kdtree_XmZp = KDTree(crd[all_edge_XmZp])

        kdtree_YpZm = KDTree(crd[all_edge_YpZm])
        kdtree_YpZp = KDTree(crd[all_edge_YpZp])
        kdtree_YmZp = KDTree(crd[all_edge_YmZp])

        offset = np.array([rve.dim_x, 0.0, 0.0])
        crd_XpYm = crd[self.edge_list_XmYm] + offset
        offset = np.array([rve.dim_x, rve.dim_y, 0.0])        
        crd_XpYp = crd[self.edge_list_XmYm] + offset   
        offset = np.array([0.0, rve.dim_y, 0.0])  
        crd_XmYp = crd[self.edge_list_XmYm] + offset       

        offset = np.array([rve.dim_x, 0.0, 0.0])
        crd_XpZm = crd[self.edge_list_XmZm] + offset
        offset = np.array([rve.dim_x, 0.0, rve.dim_z])        
        crd_XpZp = crd[self.edge_list_XmZm] + offset   
        offset = np.array([0.0, 0.0, rve.dim_z])  
        crd_XmZp = crd[self.edge_list_XmZm] + offset               

        offset = np.array([0.0, rve.dim_y, 0.0])
        crd_YpZm = crd[self.edge_list_YmZm] + offset
        offset = np.array([0.0, rve.dim_y, rve.dim_z])        
        crd_YpZp = crd[self.edge_list_YmZm] + offset   
        offset = np.array([0.0, 0.0, rve.dim_z])  
        crd_YmZp = crd[self.edge_list_YmZm] + offset

        dist_XpYm, index_XpYm = kdtree_XpYm.query(crd_XpYm, 2)
        dist_XpYp, index_XpYp = kdtree_XpYp.query(crd_XpYp, 2)
        dist_XmYp, index_XmYp = kdtree_XmYp.query(crd_XmYp, 2)

        dist_XpZm, index_XpZm = kdtree_XpZm.query(crd_XpZm, 2)
        dist_XpZp, index_XpZp = kdtree_XpZp.query(crd_XpZp, 2)
        dist_XmZp, index_XmZp = kdtree_XmZp.query(crd_XmZp, 2)       

        dist_YpZm, index_YpZm = kdtree_YpZm.query(crd_YpZm, 2)
        dist_YpZp, index_YpZp = kdtree_YpZp.query(crd_YpZp, 2)
        dist_YmZp, index_YmZp = kdtree_YmZp.query(crd_YmZp, 2)

        return (all_edge_XpYm[index_XpYm], dist_XpYm,
                all_edge_XpYp[index_XpYp], dist_XpYp,
                all_edge_XmYp[index_XmYp], dist_XmYp,
                all_edge_XpZm[index_XpZm], dist_XpZm,
                all_edge_XpZp[index_XpZp], dist_XpZp,
                all_edge_XmZp[index_XmZp], dist_XmZp,
                all_edge_YpZm[index_YpZm], dist_YpZm,
                all_edge_YpZp[index_YpZp], dist_YpZp,
                all_edge_YmZp[index_YmZp], dist_YmZp,
            )

#    def _translate_faces(self, axis, rve: Rve)
        
        
        

    # def find_neighbours(
    #     self,
    #     rve: Rve,
    #     tol: float = 1.e-8,
    # ) -> None:
        

       