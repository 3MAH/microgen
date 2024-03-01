"""
Cubic mesh for FE
"""

import warnings
from typing import Optional, NamedTuple

import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.spatial import KDTree

from .rve import Rve
from .single_mesh import SingleMesh, check_if_only_linear_tetrahedral

# try:
#     import trimesh
#     import rtree
#     import pyembree
#     USE_MULTI_RAY = True
# except ImportError:
USE_MULTI_RAY = False


class ClosestCellsOnBoundaries(NamedTuple):
    cells_id: npt.NDArray[np.int_]
    intersection_point_coords: npt.NDArray[np.float_]
    closest_opposing_cells_id: list[npt.NDArray[np.int_]]


class BoxMesh(SingleMesh):
    """
    Class to manage list of Nodes and Elements inside a Rve
    :param nodes_coords: list of nodes (np.ndarray)
    :param elements : list of elements (np.ndarray)
    :param elm_type : type of elements (np.ndarray)
    :param nodes_indices : index of node list (if different from the natural index of nodes array)
    """

    def __init__(
            self,
            nodes_coords: np.ndarray,
            elements: dict[int, npt.NDArray[np.int_]],
            pvmesh: Optional[pv.UnstructuredGrid] = None,
            nodes_indices: Optional[npt.NDArray[np.int_]] = None,
    ) -> None:
        super().__init__(
            nodes_coords=nodes_coords,
            elements=elements,
            pvmesh=pvmesh,
            nodes_indices=nodes_indices,
        )

        self._rve = None
        self._closest_points_on_boundaries = None

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

    @staticmethod
    def from_pyvista(pvmesh: pv.UnstructuredGrid):
        """Build a BoxMesh from a pyvista UnstructuredGrid or PolyData mesh (in this last case,
        the PolyData object is cast into an Unstructured Grid).
        Node and element data are not copied (by default a shallow copy is operated)
        Mesh with multiple element type are not handled.
        Note that for now singleMesh works only considering tetrahedral elements

        :param pvmesh: the mesh as a pyvista UnstructuredGrid object

        :return : A BoxMesh Object from the tetrahedral elements of the pyvista Unstructured Grid
        """
        if isinstance(pvmesh, pv.PolyData):
            pvmesh = pvmesh.cast_to_unstructured_grid()

        try:
            # extract only the tetrahedral elements
            check_if_only_linear_tetrahedral(pvmesh)
            elements = {10: pvmesh.cells_dict[10]}
            return BoxMesh(pvmesh.points, elements, pvmesh)
        except ValueError as e:
            print(e)
        return None

    def construct(
            self,
            rve: Optional[Rve] = None,
            tol: Optional[float] = 1.0e-8,
    ) -> None:
        """
        Construct a box Mesh with list of points in faces (excluding edges), edges (excluding corners) and corners.

        :param rve: RVE of the box microgen.Rve
        :param tol: tolerance to find the points on a face or an edge
        """
        if isinstance(rve, Rve) == False:
            rve = self.rve

        crd = self.nodes_coords
        closest_point_to_rve_center = np.linalg.norm(crd - rve.center, axis=1).argmin()
        self.center = crd[closest_point_to_rve_center]

        self.face_list_Xm = np.where(np.abs(crd[:, 0] - rve.x_min) < tol)[0]
        self.face_list_Xp = np.where(np.abs(crd[:, 0] - rve.x_max) < tol)[0]

        self.face_list_Ym = np.where(np.abs(crd[:, 1] - rve.y_min) < tol)[0]
        self.face_list_Yp = np.where(np.abs(crd[:, 1] - rve.y_max) < tol)[0]

        self.face_list_Zm = np.where(np.abs(crd[:, 2] - rve.z_min) < tol)[0]
        self.face_list_Zp = np.where(np.abs(crd[:, 2] - rve.z_max) < tol)[0]

        self.edge_list_XmYm = np.intersect1d(
            self.face_list_Xm, self.face_list_Ym, assume_unique=True
        )
        self.edge_list_XpYm = np.intersect1d(
            self.face_list_Xp, self.face_list_Ym, assume_unique=True
        )
        self.edge_list_XpYp = np.intersect1d(
            self.face_list_Xp, self.face_list_Yp, assume_unique=True
        )
        self.edge_list_XmYp = np.intersect1d(
            self.face_list_Xm, self.face_list_Yp, assume_unique=True
        )

        self.edge_list_XmZm = np.intersect1d(
            self.face_list_Xm, self.face_list_Zm, assume_unique=True
        )
        self.edge_list_XpZm = np.intersect1d(
            self.face_list_Xp, self.face_list_Zm, assume_unique=True
        )
        self.edge_list_XpZp = np.intersect1d(
            self.face_list_Xp, self.face_list_Zp, assume_unique=True
        )
        self.edge_list_XmZp = np.intersect1d(
            self.face_list_Xm, self.face_list_Zp, assume_unique=True
        )

        self.edge_list_YmZm = np.intersect1d(
            self.face_list_Ym, self.face_list_Zm, assume_unique=True
        )
        self.edge_list_YpZm = np.intersect1d(
            self.face_list_Yp, self.face_list_Zm, assume_unique=True
        )
        self.edge_list_YpZp = np.intersect1d(
            self.face_list_Yp, self.face_list_Zp, assume_unique=True
        )
        self.edge_list_YmZp = np.intersect1d(
            self.face_list_Ym, self.face_list_Zp, assume_unique=True
        )

        # extract corners from the intersection of edges
        self.corner_list_XmYmZm = np.intersect1d(
            self.edge_list_XmYm, self.edge_list_YmZm, assume_unique=True
        )
        self.corner_list_XmYpZm = np.intersect1d(
            self.edge_list_XmYp, self.edge_list_YpZm, assume_unique=True
        )
        self.corner_list_XpYmZm = np.intersect1d(
            self.edge_list_XpYm, self.edge_list_YmZm, assume_unique=True
        )
        self.corner_list_XpYpZm = np.intersect1d(
            self.edge_list_XpYp, self.edge_list_YpZm, assume_unique=True
        )
        self.corner_list_XmYmZp = np.intersect1d(
            self.edge_list_XmYm, self.edge_list_YmZp, assume_unique=True
        )
        self.corner_list_XmYpZp = np.intersect1d(
            self.edge_list_XmYp, self.edge_list_YpZp, assume_unique=True
        )
        self.corner_list_XpYmZp = np.intersect1d(
            self.edge_list_XpYm, self.edge_list_YmZp, assume_unique=True
        )
        self.corner_list_XpYpZp = np.intersect1d(
            self.edge_list_XpYp, self.edge_list_YpZp, assume_unique=True
        )

        # Remove nodes that belong to several sets
        all_corners = np.hstack(
            (
                self.corner_list_XmYmZm,
                self.corner_list_XmYpZm,
                self.corner_list_XpYmZm,
                self.corner_list_XpYpZm,
                self.corner_list_XmYmZp,
                self.corner_list_XmYpZp,
                self.corner_list_XpYmZp,
                self.corner_list_XpYpZp,
            )
        )

        self.edge_list_XmYm = np.setdiff1d(
            self.edge_list_XmYm, all_corners, assume_unique=True
        )
        self.edge_list_XpYm = np.setdiff1d(
            self.edge_list_XpYm, all_corners, assume_unique=True
        )
        self.edge_list_XpYp = np.setdiff1d(
            self.edge_list_XpYp, all_corners, assume_unique=True
        )
        self.edge_list_XmYp = np.setdiff1d(
            self.edge_list_XmYp, all_corners, assume_unique=True
        )

        self.edge_list_XmZm = np.setdiff1d(
            self.edge_list_XmZm, all_corners, assume_unique=True
        )
        self.edge_list_XpZm = np.setdiff1d(
            self.edge_list_XpZm, all_corners, assume_unique=True
        )
        self.edge_list_XpZp = np.setdiff1d(
            self.edge_list_XpZp, all_corners, assume_unique=True
        )
        self.edge_list_XmZp = np.setdiff1d(
            self.edge_list_XmZp, all_corners, assume_unique=True
        )

        self.edge_list_YmZm = np.setdiff1d(
            self.edge_list_YmZm, all_corners, assume_unique=True
        )
        self.edge_list_YpZm = np.setdiff1d(
            self.edge_list_YpZm, all_corners, assume_unique=True
        )
        self.edge_list_YpZp = np.setdiff1d(
            self.edge_list_YpZp, all_corners, assume_unique=True
        )
        self.edge_list_YmZp = np.setdiff1d(
            self.edge_list_YmZp, all_corners, assume_unique=True
        )

        all_edges_corners = np.hstack(
            (
                self.edge_list_XmYm,
                self.edge_list_XpYm,
                self.edge_list_XpYp,
                self.edge_list_XmYp,
                self.edge_list_XmZm,
                self.edge_list_XpZm,
                self.edge_list_XpZp,
                self.edge_list_XmZp,
                self.edge_list_YmZm,
                self.edge_list_YpZm,
                self.edge_list_YpZp,
                self.edge_list_YmZp,
                all_corners,
            )
        )

        self.face_list_Xm = np.setdiff1d(
            self.face_list_Xm, all_edges_corners, assume_unique=True
        )
        self.face_list_Xp = np.setdiff1d(
            self.face_list_Xp, all_edges_corners, assume_unique=True
        )
        self.face_list_Ym = np.setdiff1d(
            self.face_list_Ym, all_edges_corners, assume_unique=True
        )
        self.face_list_Yp = np.setdiff1d(
            self.face_list_Yp, all_edges_corners, assume_unique=True
        )
        self.face_list_Zm = np.setdiff1d(
            self.face_list_Zm, all_edges_corners, assume_unique=True
        )
        self.face_list_Zp = np.setdiff1d(
            self.face_list_Zp, all_edges_corners, assume_unique=True
        )

        self.edge_list_XmYm = self.edge_list_XmYm[
            np.argsort(crd[self.edge_list_XmYm, 2])
        ]
        self.edge_list_XpYm = self.edge_list_XpYm[
            np.argsort(crd[self.edge_list_XpYm, 2])
        ]
        self.edge_list_XpYp = self.edge_list_XpYp[
            np.argsort(crd[self.edge_list_XpYp, 2])
        ]
        self.edge_list_XmYp = self.edge_list_XmYp[
            np.argsort(crd[self.edge_list_XmYp, 2])
        ]

        self.edge_list_XmZm = self.edge_list_XmZm[
            np.argsort(crd[self.edge_list_XmZm, 1])
        ]
        self.edge_list_XpZm = self.edge_list_XpZm[
            np.argsort(crd[self.edge_list_XpZm, 1])
        ]
        self.edge_list_XpZp = self.edge_list_XpZp[
            np.argsort(crd[self.edge_list_XpZp, 1])
        ]
        self.edge_list_XmZp = self.edge_list_XmZp[
            np.argsort(crd[self.edge_list_XmZp, 1])
        ]

        self.edge_list_YmZm = self.edge_list_YmZm[
            np.argsort(crd[self.edge_list_YmZm, 0])
        ]
        self.edge_list_YpZm = self.edge_list_YpZm[
            np.argsort(crd[self.edge_list_YpZm, 0])
        ]
        self.edge_list_YpZp = self.edge_list_YpZp[
            np.argsort(crd[self.edge_list_YpZp, 0])
        ]
        self.edge_list_YmZp = self.edge_list_YmZp[
            np.argsort(crd[self.edge_list_YmZp, 0])
        ]

        decimal_round = int(-np.log10(tol) - 1)
        self.face_list_Xm = self.face_list_Xm[
            np.lexsort(
                (
                    crd[self.face_list_Xm, 1],
                    crd[self.face_list_Xm, 2].round(decimal_round),
                )
            )
        ]
        self.face_list_Xp = self.face_list_Xp[
            np.lexsort(
                (
                    crd[self.face_list_Xp, 1],
                    crd[self.face_list_Xp, 2].round(decimal_round),
                )
            )
        ]
        self.face_list_Ym = self.face_list_Ym[
            np.lexsort(
                (
                    crd[self.face_list_Ym, 0],
                    crd[self.face_list_Ym, 2].round(decimal_round),
                )
            )
        ]
        self.face_list_Yp = self.face_list_Yp[
            np.lexsort(
                (
                    crd[self.face_list_Yp, 0],
                    crd[self.face_list_Yp, 2].round(decimal_round),
                )
            )
        ]
        self.face_list_Zm = self.face_list_Zm[
            np.lexsort(
                (
                    crd[self.face_list_Zm, 0],
                    crd[self.face_list_Zm, 1].round(decimal_round),
                )
            )
        ]
        self.face_list_Zp = self.face_list_Zp[
            np.lexsort(
                (
                    crd[self.face_list_Zp, 0],
                    crd[self.face_list_Zp, 1].round(decimal_round),
                )
            )
        ]

        self.corners = [
            self.corner_list_XmYmZm,
            self.corner_list_XmYpZm,
            self.corner_list_XpYmZm,
            self.corner_list_XpYpZm,
            self.corner_list_XmYmZp,
            self.corner_list_XmYpZp,
            self.corner_list_XpYmZp,
            self.corner_list_XpYpZp,
        ]

        self.edges = [
            self.edge_list_XmYm,
            self.edge_list_XpYm,
            self.edge_list_XpYp,
            self.edge_list_XmYp,
            self.edge_list_XmZm,
            self.edge_list_XpZm,
            self.edge_list_XpZp,
            self.edge_list_XmZp,
            self.edge_list_YmZm,
            self.edge_list_YpZm,
            self.edge_list_YpZp,
            self.edge_list_YmZp,
        ]

        self.faces = [
            self.face_list_Xm,
            self.face_list_Xp,
            self.face_list_Ym,
            self.face_list_Yp,
            self.face_list_Zm,
            self.face_list_Zp,
        ]

    @property
    def rve(
            self,
    ) -> Rve:
        """
        Return a representative volume element (Rve) of the considered mesh

        :return Rve: a Rve object attached to the considered BoxMesh object
        """

        if isinstance(self._rve, Rve):
            return self._rve
        else:
            return self._build_rve()

    @rve.setter
    def rve(
            self,
            value: Rve,
    ) -> None:
        """
        set a representative volume element (Rve) of the considered mesh

        :param value: Rve of the box to set
        """

        self._rve = value

    def _build_rve(
            self,
    ) -> Rve:
        """
        build a representative volume element (Rve) of the considered mesh from its bounding box

        :return rve: RVE of the mesh bounding box
        """
        if not (isinstance(self._pvmesh, pv.UnstructuredGrid)):
            self._mesh = self.to_pyvista()
        xmin, xmax, ymin, ymax, zmin, zmax = self._mesh.bounds
        return Rve.from_min_max(xmin, xmax, ymin, ymax, zmin, zmax)

    def _closest_points_on_faces(
            self,
            k_neighbours: int = 3,
            rve: Optional[Rve] = None,
            tol: Optional[float] = 1.0e-8,
    ) -> dict[str, tuple[npt.NDArray[np.int_], npt.NDArray[np.float_]]]:
        """
        Find the closest points on opposite face to write interpolation relationship
        if a displacement condition between pair nodes is defined on such opposite surfaces
        It takes the set of points on faces 'face_Xm', 'face_Ym', 'face_Zm' (excluding edges and corners)
        and returns for each opposite face the closest points (index of the points, distance) of the
        projected point on the corresponding opposite face

        :param k_neighbours : number of closest points
        :param rve : RVE of the mesh bounding box. if None, the rve is built from the mesh bounding box

        :return dict: a dictionnary with (np.array) of indices and a np.array of distances for each neighbor:
            'face_Xp' : (index[0], dist[0]),
            'face_Yp' : (index[1], dist[1]),
            'face_Zp' : (index[2], dist[2])
        """

        if isinstance(rve, Rve) == False:
            rve = self.rve

        crd = self.nodes_coords

        all_face_Xp = np.hstack(
            (
                self.face_list_Xp,
                self.edge_list_XpYm,
                self.edge_list_XpYp,
                self.edge_list_XpZm,
                self.edge_list_XpZp,
                self.corner_list_XpYmZm,
                self.corner_list_XpYpZm,
                self.corner_list_XpYmZp,
                self.corner_list_XpYpZp,
            )
        )
        all_face_Yp = np.hstack(
            (
                self.face_list_Yp,
                self.edge_list_YpZm,
                self.edge_list_YpZp,
                self.edge_list_XmYp,
                self.edge_list_XpYp,
                self.corner_list_XmYpZm,
                self.corner_list_XpYpZm,
                self.corner_list_XmYpZp,
                self.corner_list_XpYpZp,
            )
        )
        all_face_Zp = np.hstack(
            (
                self.face_list_Zp,
                self.edge_list_YmZp,
                self.edge_list_YpZp,
                self.edge_list_XmZp,
                self.edge_list_XpZp,
                self.corner_list_XmYmZp,
                self.corner_list_XpYmZp,
                self.corner_list_XmYpZp,
                self.corner_list_XpYpZp,
            )
        )

        kdTrees = [
            KDTree(crd[all_face_Xp]),
            KDTree(crd[all_face_Yp]),
            KDTree(crd[all_face_Zp]),
        ]

        offsets = [
            np.array([rve.dim_x, 0.0, 0.0]),
            np.array([0.0, rve.dim_y, 0.0]),
            np.array([0.0, 0.0, rve.dim_z]),
        ]

        faces_m = [self.face_list_Xm, self.face_list_Ym, self.face_list_Zm]
        all_faces_p = [all_face_Xp, all_face_Yp, all_face_Zp]

        dist = []
        index = []

        for i, face in enumerate(faces_m):
            offset = offsets[i]
            crd_face = crd[face] + offset
            minimum_query_points = min(len(all_faces_p[i]), k_neighbours)

            if minimum_query_points < k_neighbours:
                warnings.warn(
                    "Number of query points is greater than the number of points in the KDTree."
                )

            dist_temp, index_temp = kdTrees[i].query(crd_face, minimum_query_points)
            all_faces_p_i = all_faces_p[i]
            index_temp_list = all_faces_p_i[index_temp].tolist()
            if k_neighbours == 1:
                dist_temp_list = [[d] for d in dist_temp]
                index_temp_list = [[i] for i in index_temp_list]
            else:
                dist_temp_list = dist_temp.tolist()

            # If pair nodes exist (opposite node exactly match), return only the pair node (the closest neighbour)
            if k_neighbours > 1:
                for i in range(0, len(dist_temp)):
                    if dist_temp_list[i][0] < tol:
                        dist_temp_list[i] = dist_temp_list[i][0:1]
                        index_temp_list[i] = index_temp_list[i][0:1]

            dist.append(dist_temp_list)
            index.append(index_temp_list)

        return {
            "face_Xp": (index[0], dist[0]),
            "face_Yp": (index[1], dist[1]),
            "face_Zp": (index[2], dist[2]),
        }

    def _closest_points_on_edges(
            self,
            rve: Rve = None,
            tol: Optional[float] = 1.0e-8,
    ) -> dict[str, tuple[npt.NDArray[np.int_], npt.NDArray[np.float_]]]:
        """
        Find the closest points on opposite edges to write interpolation relationship
        if a displacement condition between pair nodes is defined on such opposite surfaces
        Note : the closest neighbours is set as 2 for edges

        :param rve : RVE of the mesh bounding box. if None, the rve is built from the mesh bounding box
        """

        if isinstance(rve, Rve) == False:
            rve = self.rve

        crd = self.nodes_coords

        all_edge_XpYm = np.hstack(
            (self.edge_list_XpYm, self.corner_list_XpYmZm, self.corner_list_XpYmZp)
        )
        all_edge_XpYp = np.hstack(
            (self.edge_list_XpYp, self.corner_list_XpYpZm, self.corner_list_XpYpZp)
        )
        all_edge_XmYp = np.hstack(
            (self.edge_list_XmYp, self.corner_list_XmYpZm, self.corner_list_XmYpZp)
        )

        all_edge_XpZm = np.hstack(
            (self.edge_list_XpZm, self.corner_list_XpYmZm, self.corner_list_XpYpZm)
        )
        all_edge_XpZp = np.hstack(
            (self.edge_list_XpZp, self.corner_list_XpYmZp, self.corner_list_XpYpZp)
        )
        all_edge_XmZp = np.hstack(
            (self.edge_list_XmZp, self.corner_list_XmYmZp, self.corner_list_XmYpZp)
        )

        all_edge_YpZm = np.hstack(
            (self.edge_list_YmZm, self.corner_list_XmYpZm, self.corner_list_XpYpZm)
        )
        all_edge_YpZp = np.hstack(
            (self.edge_list_YpZp, self.corner_list_XmYpZp, self.corner_list_XpYpZp)
        )
        all_edge_YmZp = np.hstack(
            (self.edge_list_YmZp, self.corner_list_XmYmZp, self.corner_list_XpYmZp)
        )

        kdTrees = [
            KDTree(crd[all_edge_XpYm]),
            KDTree(crd[all_edge_XpYp]),
            KDTree(crd[all_edge_XmYp]),
            KDTree(crd[all_edge_XpZm]),
            KDTree(crd[all_edge_XpZp]),
            KDTree(crd[all_edge_XmZp]),
            KDTree(crd[all_edge_YpZm]),
            KDTree(crd[all_edge_YpZp]),
            KDTree(crd[all_edge_YmZp]),
        ]

        offsets = [
            np.array([rve.dim_x, 0.0, 0.0]),
            np.array([rve.dim_x, rve.dim_y, 0.0]),
            np.array([0.0, rve.dim_y, 0.0]),
            np.array([rve.dim_x, 0.0, 0.0]),
            np.array([rve.dim_x, 0.0, rve.dim_z]),
            np.array([0.0, 0.0, rve.dim_z]),
            np.array([0.0, rve.dim_y, 0.0]),
            np.array([0.0, rve.dim_y, rve.dim_z]),
            np.array([0.0, 0.0, rve.dim_z]),
        ]

        edges_m = [
            self.edge_list_XmYm,
            self.edge_list_XmYm,
            self.edge_list_XmYm,
            self.edge_list_XmZm,
            self.edge_list_XmZm,
            self.edge_list_XmZm,
            self.edge_list_YmZm,
            self.edge_list_YmZm,
            self.edge_list_YmZm,
        ]

        all_edges_p = [
            all_edge_XpYm,
            all_edge_XpYp,
            all_edge_XmYp,
            all_edge_XpZm,
            all_edge_XpZp,
            all_edge_XmZp,
            all_edge_YpZm,
            all_edge_YpZp,
            all_edge_YmZp,
        ]

        dist = []
        index = []

        for i, edge in enumerate(edges_m):
            offset = offsets[i]
            crd_edge = crd[edge] + offset

            dist_temp, index_temp = kdTrees[i].query(crd_edge, 2)

            dist_temp_list = dist_temp.tolist()
            all_edges_p_i = all_edges_p[i]
            index_temp_list = all_edges_p_i[index_temp].tolist()

            for i in range(0, len(dist_temp)):
                if dist_temp_list[i][0] < tol:
                    dist_temp_list[i] = dist_temp_list[i][0:1]
                    index_temp_list[i] = index_temp_list[i][0:1]

            dist.append(dist_temp_list)
            index.append(index_temp_list)

        return {
            "edge_XpYm": (index[0], dist[0]),
            "edge_XpYp": (index[1], dist[1]),
            "edge_XmYp": (index[2], dist[2]),
            "edge_XpZm": (index[3], dist[3]),
            "edge_XpZp": (index[4], dist[4]),
            "edge_XmZp": (index[5], dist[5]),
            "edge_YpZm": (index[6], dist[6]),
            "edge_YpZp": (index[7], dist[7]),
            "edge_YmZp": (index[8], dist[8]),
        }

    def closest_points_on_boundaries(
            self,
            k_neighbours: int = 3,
            rve: Optional[Rve] = None,
            tol: Optional[float] = 1.0e-8,
    ) -> dict[str, tuple[npt.NDArray[np.int_], npt.NDArray[np.float_]]]:
        """
        Find the closest points on faces and edges to write interpolation relationship
        if a displacement condition between pair nodes is defined

        :param rve : RVE of the mesh bounding box. if None, the rve is built from the mesh bounding box
        """

        if isinstance(self._closest_points_on_boundaries, dict):
            return self._closest_points_on_boundaries
        else:
            if isinstance(rve, Rve) == False:
                rve = self.rve

            dict_faces = self._closest_points_on_faces(k_neighbours, rve, tol)
            dict_edges = self._closest_points_on_edges(rve, tol)

            self._dict_faces = dict_faces
            self._dict_edges = dict_edges

            self._closest_points_on_boundaries = {**dict_faces, **dict_edges}
            return self._closest_points_on_boundaries

    def boundary_elements(
            self,
            rve: Optional[Rve] = None,
            tol: Optional[float] = 1.0e-4,
    ) -> tuple[pv.PolyData, npt.NDArray[np.int_]]:
        if isinstance(rve, Rve) == False:
            rve = self.rve

        normals = [
            np.array([-1.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, -1.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, -1.0]),
            np.array([0.0, 0.0, 1.0]),
        ]
        origins_p = [
            np.array([rve.x_min, rve.center[1], rve.center[2]]),
            np.array([rve.x_max, rve.center[1], rve.center[2]]),
            np.array([rve.center[0], rve.y_min, rve.center[2]]),
            np.array([rve.center[0], rve.y_max, rve.center[2]]),
            np.array([rve.center[0], rve.center[1], rve.z_min]),
            np.array([rve.center[0], rve.center[1], rve.z_max]),
        ]
        size_planes = [
            2.0 * rve.dx,
            2.0 * rve.dx,
            2.0 * rve.dy,
            2.0 * rve.dy,
            2.0 * rve.dz,
            2.0 * rve.dz,
        ]

        surface = self.surface
        surface["CellIDs"] = np.arange(surface.n_cells)

        boundary_elements = pv.PolyData()

        for normal, origin_p, size_plane in zip(normals, origins_p, size_planes):
            plane = pv.Plane(
                center=origin_p, direction=normal, i_size=size_plane, j_size=size_plane
            )
            surface.compute_implicit_distance(plane, inplace=True)
            surface_p = surface.threshold(
                [-tol, tol], all_scalars=True, scalars="implicit_distance"
            ).extract_surface()
            boundary_elements = boundary_elements.append_polydata(surface_p)

        return boundary_elements, boundary_elements["CellIDs"]

    def closest_cells_on_boundaries(
            self,
            rve: Optional[Rve] = None,
            tol: Optional[float] = 1.0e-8,
    ) -> dict[str, ClosestCellsOnBoundaries]:
        """
        Find the cells to which a given point belongs to while using a ray tracing normal to the face on which it belongs
        :param rve : RVE of the mesh bounding box. if None, the rve is built from the mesh bounding box
        :param tol : tolerance to evaluate the threshold between the cells and the boundary RVE of the mesh bounding box. if None, the rve is built from the mesh bounding box
        """

        if isinstance(rve, Rve) == False:
            rve = self.rve

        crd = self.nodes_coords

        normals = [
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        ]
        origins_p = [
            np.array([rve.x_max, rve.center[1], rve.center[2]]),
            np.array([rve.center[0], rve.y_max, rve.center[2]]),
            np.array([rve.center[0], rve.center[1], rve.z_max]),
        ]
        size_planes = [2.0 * rve.dx, 2.0 * rve.dy, 2.0 * rve.dz]

        faces_m = [
            crd[self.face_list_Xm],
            crd[self.face_list_Ym],
            crd[self.face_list_Zm],
        ]
        surface = self.surface
        surface["CellIDs"] = np.arange(surface.n_cells)

        list_cells_for_each_face = []
        list_ray_trace = []

        for i, face in enumerate(faces_m):
            plane = pv.Plane(
                center=origins_p[i],
                direction=normals[i],
                i_size=size_planes[i],
                j_size=size_planes[i],
            )
            surface.compute_implicit_distance(plane, inplace=True)

            surface_p = surface.threshold(
                [-tol, tol], all_scalars=True, scalars="implicit_distance"
            ).extract_surface()
            list_cells_for_each_face.append(surface_p["CellIDs"])

            directions = np.tile(normals[i], (np.shape(faces_m[i])[0], 1))

            raytraceresult = []
            if USE_MULTI_RAY:
                raytraceresult = surface_p.multi_ray_trace(
                    origins=faces_m[i], directions=directions
                )
            else:
                intersection_points = []
                intersection_rays = []
                intersection_cells = []
                for j, face_m_i in enumerate(faces_m[i]):
                    end_point = face_m_i + size_planes[i] * directions[j]
                    intersection_point, intersection_cell = surface_p.ray_trace(
                        origin=face_m_i, end_point=end_point
                    )
                    intersection_ray = np.full_like(intersection_cell, j)
                    intersection_points.extend(intersection_point.tolist())
                    intersection_rays.extend(intersection_ray.tolist())
                    intersection_cells.extend(intersection_cell.tolist())
                raytraceresult = (
                    np.asarray(intersection_points),
                    np.asarray(intersection_rays),
                    np.asarray(intersection_cells),
                )

            list_ray_trace.append(raytraceresult)

        return {
            "face_Xp": ClosestCellsOnBoundaries(
                list_cells_for_each_face[0],
                list_ray_trace[0][0],
                list_cells_for_each_face[0][list_ray_trace[0][2]],
            ),
            "face_Yp": ClosestCellsOnBoundaries(
                list_cells_for_each_face[1],
                list_ray_trace[1][0],
                list_cells_for_each_face[1][list_ray_trace[1][2]],
            ),
            "face_Zp": ClosestCellsOnBoundaries(
                list_cells_for_each_face[2],
                list_ray_trace[2][0],
                list_cells_for_each_face[2][list_ray_trace[2][2]],
            ),
        }
