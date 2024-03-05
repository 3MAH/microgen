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

# We could in future versions make benefit of the embree library
# with multi ray tracing embedded in PyVista.
# However, a very old version of embree (2) is nowadays bound in python
# try:
#     import trimesh
#     import rtree
#     import pyembree
#     USE_MULTI_RAY = True
# except ImportError:

USE_MULTI_RAY = False


class ClosestCellsOnBoundaries(NamedTuple):
    """
    Class to manage closest cells on boundaries
    :param cells_id: list of cells identification number
    :param intersection_point_coords: list of intersection points of projected
    rays from one face's nodes with the opposing face
    :param closest_opposing_cells_id: list of closest cells identification number on opposing face
    """
    cells_id: npt.NDArray[np.int_]
    intersection_point_coords: npt.NDArray[np.float_]
    closest_opposing_cells_id: list[npt.NDArray[np.int_]]


class BoxMesh(SingleMesh):
    """
    Class to manage list of Nodes and Elements inside a Rve
    :param nodes_coords: list of nodes (np.ndarray)
    :param elements: list of elements (np.ndarray)
    :param nodes_indices : index of node list (if different from the natural index of nodes array)
    """

    def __init__(
            self,
            nodes_coords: npt.NDArray[np.float_],
            elements: dict[pv.CellType, npt.NDArray[np.int_]],
            pvmesh: Optional[pv.UnstructuredGrid] = None,
            nodes_indices: Optional[npt.NDArray[np.int_]] = None,
    ) -> None:
        super().__init__(
            nodes_coords=nodes_coords,
            elements=elements,
            pvmesh=pvmesh,
            nodes_indices=nodes_indices,
        )

        self._rve: Optional[Rve] = None
        self._closest_points_on_boundaries: Optional[
            dict[str, tuple[npt.NDArray[np.int_], npt.NDArray[np.float_]]]] = None

        self.center: Optional[npt.NDArray[np.float_]] = None

        self._construct(self.rve)

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

        # extract only the tetrahedral elements
        # raises an Exception if 3d elements other than linear tetrahedra are found in the mesh
        check_if_only_linear_tetrahedral(pvmesh)
        elements = {pv.CellType.TETRA: pvmesh.cells_dict[pv.CellType.TETRA]}
        return BoxMesh(pvmesh.points, elements, pvmesh)

    def _construct(
            self,
            rve: Optional[Rve] = None,
            tol: float = 1.0e-8,
    ) -> None:
        """
        Construct a box Mesh with list of points in faces (excluding edges), edges (excluding corners) and corners.

        :param rve: RVE of the box microgen.Rve
        :param tol: tolerance to find the points on a face or an edge
        """
        if rve is None:
            rve = self.rve

        crd = self.nodes_coords
        closest_point_to_rve_center = np.linalg.norm(crd - rve.center, axis=1).argmin()
        self.center = crd[closest_point_to_rve_center]

        self.face_list_xm = np.where(np.abs(crd[:, 0] - rve.x_min) < tol)[0]
        self.face_list_xp = np.where(np.abs(crd[:, 0] - rve.x_max) < tol)[0]

        self.face_list_ym = np.where(np.abs(crd[:, 1] - rve.y_min) < tol)[0]
        self.face_list_yp = np.where(np.abs(crd[:, 1] - rve.y_max) < tol)[0]

        self.face_list_zm = np.where(np.abs(crd[:, 2] - rve.z_min) < tol)[0]
        self.face_list_zp = np.where(np.abs(crd[:, 2] - rve.z_max) < tol)[0]

        self.edge_list_xm_ym = np.intersect1d(
            self.face_list_xm, self.face_list_ym, assume_unique=True
        )
        self.edge_list_xp_ym = np.intersect1d(
            self.face_list_xp, self.face_list_ym, assume_unique=True
        )
        self.edge_list_xp_yp = np.intersect1d(
            self.face_list_xp, self.face_list_yp, assume_unique=True
        )
        self.edge_list_xm_yp = np.intersect1d(
            self.face_list_xm, self.face_list_yp, assume_unique=True
        )

        self.edge_list_xm_zm = np.intersect1d(
            self.face_list_xm, self.face_list_zm, assume_unique=True
        )
        self.edge_list_xp_zm = np.intersect1d(
            self.face_list_xp, self.face_list_zm, assume_unique=True
        )
        self.edge_list_xp_zp = np.intersect1d(
            self.face_list_xp, self.face_list_zp, assume_unique=True
        )
        self.edge_list_xm_zp = np.intersect1d(
            self.face_list_xm, self.face_list_zp, assume_unique=True
        )

        self.edge_list_ym_zm = np.intersect1d(
            self.face_list_ym, self.face_list_zm, assume_unique=True
        )
        self.edge_list_yp_zm = np.intersect1d(
            self.face_list_yp, self.face_list_zm, assume_unique=True
        )
        self.edge_list_yp_zp = np.intersect1d(
            self.face_list_yp, self.face_list_zp, assume_unique=True
        )
        self.edge_list_ym_zp = np.intersect1d(
            self.face_list_ym, self.face_list_zp, assume_unique=True
        )

        # extract corners from the intersection of edges
        self.corner_list_xm_ym_zm = np.intersect1d(
            self.edge_list_xm_ym, self.edge_list_ym_zm, assume_unique=True
        )
        self.corner_list_xm_yp_zm = np.intersect1d(
            self.edge_list_xm_yp, self.edge_list_yp_zm, assume_unique=True
        )
        self.corner_list_xp_ym_zm = np.intersect1d(
            self.edge_list_xp_ym, self.edge_list_ym_zm, assume_unique=True
        )
        self.corner_list_xp_yp_zm = np.intersect1d(
            self.edge_list_xp_yp, self.edge_list_yp_zm, assume_unique=True
        )
        self.corner_list_xm_ym_zp = np.intersect1d(
            self.edge_list_xm_ym, self.edge_list_ym_zp, assume_unique=True
        )
        self.corner_list_xm_yp_zp = np.intersect1d(
            self.edge_list_xm_yp, self.edge_list_yp_zp, assume_unique=True
        )
        self.corner_list_xp_ym_zp = np.intersect1d(
            self.edge_list_xp_ym, self.edge_list_ym_zp, assume_unique=True
        )
        self.corner_list_xp_yp_zp = np.intersect1d(
            self.edge_list_xp_yp, self.edge_list_yp_zp, assume_unique=True
        )

        # Remove nodes that belong to several sets
        all_corners = np.hstack(
            (
                self.corner_list_xm_ym_zm,
                self.corner_list_xm_yp_zm,
                self.corner_list_xp_ym_zm,
                self.corner_list_xp_yp_zm,
                self.corner_list_xm_ym_zp,
                self.corner_list_xm_yp_zp,
                self.corner_list_xp_ym_zp,
                self.corner_list_xp_yp_zp,
            )
        )

        self.edge_list_xm_ym = np.setdiff1d(
            self.edge_list_xm_ym, all_corners, assume_unique=True
        )
        self.edge_list_xp_ym = np.setdiff1d(
            self.edge_list_xp_ym, all_corners, assume_unique=True
        )
        self.edge_list_xp_yp = np.setdiff1d(
            self.edge_list_xp_yp, all_corners, assume_unique=True
        )
        self.edge_list_xm_yp = np.setdiff1d(
            self.edge_list_xm_yp, all_corners, assume_unique=True
        )

        self.edge_list_xm_zm = np.setdiff1d(
            self.edge_list_xm_zm, all_corners, assume_unique=True
        )
        self.edge_list_xp_zm = np.setdiff1d(
            self.edge_list_xp_zm, all_corners, assume_unique=True
        )
        self.edge_list_xp_zp = np.setdiff1d(
            self.edge_list_xp_zp, all_corners, assume_unique=True
        )
        self.edge_list_xm_zp = np.setdiff1d(
            self.edge_list_xm_zp, all_corners, assume_unique=True
        )

        self.edge_list_ym_zm = np.setdiff1d(
            self.edge_list_ym_zm, all_corners, assume_unique=True
        )
        self.edge_list_yp_zm = np.setdiff1d(
            self.edge_list_yp_zm, all_corners, assume_unique=True
        )
        self.edge_list_yp_zp = np.setdiff1d(
            self.edge_list_yp_zp, all_corners, assume_unique=True
        )
        self.edge_list_ym_zp = np.setdiff1d(
            self.edge_list_ym_zp, all_corners, assume_unique=True
        )

        all_edges_corners = np.hstack(
            (
                self.edge_list_xm_ym,
                self.edge_list_xp_ym,
                self.edge_list_xp_yp,
                self.edge_list_xm_yp,
                self.edge_list_xm_zm,
                self.edge_list_xp_zm,
                self.edge_list_xp_zp,
                self.edge_list_xm_zp,
                self.edge_list_ym_zm,
                self.edge_list_yp_zm,
                self.edge_list_yp_zp,
                self.edge_list_ym_zp,
                all_corners,
            )
        )

        self.face_list_xm = np.setdiff1d(
            self.face_list_xm, all_edges_corners, assume_unique=True
        )
        self.face_list_xp = np.setdiff1d(
            self.face_list_xp, all_edges_corners, assume_unique=True
        )
        self.face_list_ym = np.setdiff1d(
            self.face_list_ym, all_edges_corners, assume_unique=True
        )
        self.face_list_yp = np.setdiff1d(
            self.face_list_yp, all_edges_corners, assume_unique=True
        )
        self.face_list_zm = np.setdiff1d(
            self.face_list_zm, all_edges_corners, assume_unique=True
        )
        self.face_list_zp = np.setdiff1d(
            self.face_list_zp, all_edges_corners, assume_unique=True
        )

        self.edge_list_xm_ym = self.edge_list_xm_ym[
            np.argsort(crd[self.edge_list_xm_ym, 2])
        ]
        self.edge_list_xp_ym = self.edge_list_xp_ym[
            np.argsort(crd[self.edge_list_xp_ym, 2])
        ]
        self.edge_list_xp_yp = self.edge_list_xp_yp[
            np.argsort(crd[self.edge_list_xp_yp, 2])
        ]
        self.edge_list_xm_yp = self.edge_list_xm_yp[
            np.argsort(crd[self.edge_list_xm_yp, 2])
        ]

        self.edge_list_xm_zm = self.edge_list_xm_zm[
            np.argsort(crd[self.edge_list_xm_zm, 1])
        ]
        self.edge_list_xp_zm = self.edge_list_xp_zm[
            np.argsort(crd[self.edge_list_xp_zm, 1])
        ]
        self.edge_list_xp_zp = self.edge_list_xp_zp[
            np.argsort(crd[self.edge_list_xp_zp, 1])
        ]
        self.edge_list_xm_zp = self.edge_list_xm_zp[
            np.argsort(crd[self.edge_list_xm_zp, 1])
        ]

        self.edge_list_ym_zm = self.edge_list_ym_zm[
            np.argsort(crd[self.edge_list_ym_zm, 0])
        ]
        self.edge_list_yp_zm = self.edge_list_yp_zm[
            np.argsort(crd[self.edge_list_yp_zm, 0])
        ]
        self.edge_list_yp_zp = self.edge_list_yp_zp[
            np.argsort(crd[self.edge_list_yp_zp, 0])
        ]
        self.edge_list_ym_zp = self.edge_list_ym_zp[
            np.argsort(crd[self.edge_list_ym_zp, 0])
        ]

        decimal_round = int(-np.log10(tol) - 1)
        self.face_list_xm = self.face_list_xm[
            np.lexsort(
                (
                    crd[self.face_list_xm, 1],
                    crd[self.face_list_xm, 2].round(decimal_round),
                )
            )
        ]
        self.face_list_xp = self.face_list_xp[
            np.lexsort(
                (
                    crd[self.face_list_xp, 1],
                    crd[self.face_list_xp, 2].round(decimal_round),
                )
            )
        ]
        self.face_list_ym = self.face_list_ym[
            np.lexsort(
                (
                    crd[self.face_list_ym, 0],
                    crd[self.face_list_ym, 2].round(decimal_round),
                )
            )
        ]
        self.face_list_yp = self.face_list_yp[
            np.lexsort(
                (
                    crd[self.face_list_yp, 0],
                    crd[self.face_list_yp, 2].round(decimal_round),
                )
            )
        ]
        self.face_list_zm = self.face_list_zm[
            np.lexsort(
                (
                    crd[self.face_list_zm, 0],
                    crd[self.face_list_zm, 1].round(decimal_round),
                )
            )
        ]
        self.face_list_zp = self.face_list_zp[
            np.lexsort(
                (
                    crd[self.face_list_zp, 0],
                    crd[self.face_list_zp, 1].round(decimal_round),
                )
            )
        ]

        self.corners = [
            self.corner_list_xm_ym_zm,
            self.corner_list_xm_yp_zm,
            self.corner_list_xp_ym_zm,
            self.corner_list_xp_yp_zm,
            self.corner_list_xm_ym_zp,
            self.corner_list_xm_yp_zp,
            self.corner_list_xp_ym_zp,
            self.corner_list_xp_yp_zp,
        ]

        self.edges = [
            self.edge_list_xm_ym,
            self.edge_list_xp_ym,
            self.edge_list_xp_yp,
            self.edge_list_xm_yp,
            self.edge_list_xm_zm,
            self.edge_list_xp_zm,
            self.edge_list_xp_zp,
            self.edge_list_xm_zp,
            self.edge_list_ym_zm,
            self.edge_list_yp_zm,
            self.edge_list_yp_zp,
            self.edge_list_ym_zp,
        ]

        self.faces = [
            self.face_list_xm,
            self.face_list_xp,
            self.face_list_ym,
            self.face_list_yp,
            self.face_list_zm,
            self.face_list_zp,
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
        if not isinstance(self._pvmesh, pv.UnstructuredGrid):
            self._mesh = self.to_pyvista()
        xmin, xmax, ymin, ymax, zmin, zmax = self._mesh.bounds
        return Rve.from_min_max(float(xmin), float(xmax), float(ymin),
                                float(ymax), float(zmin), float(zmax))

    def _closest_points_on_faces(
            self,
            k_neighbours: int = 3,
            rve: Optional[Rve] = None,
            tol: float = 1.0e-8,
    ) -> dict[str, tuple[npt.NDArray[np.int_], npt.NDArray[np.float_]]]:
        """
        Find the closest points on opposite face to write interpolation relationship
        if a displacement condition between pair nodes is defined on such opposite surfaces
        It takes the set of points on faces 'face_Xm', 'face_Ym', 'face_Zm' (excluding edges and corners)
        and returns for each opposite face the closest points (index of the points, distance) of the
        projected point on the corresponding opposite face

        :param k_neighbours : number of closest points
        :param rve : RVE of the mesh bounding box. if None, the rve is built from the mesh bounding box
        :param tol: tolerance

        :return dict: a dictionary with (np.array) of indices and a np.array of distances for each neighbor:
            'face_Xp' : (index[0], dist[0]),
            'face_Yp' : (index[1], dist[1]),
            'face_Zp' : (index[2], dist[2])
        """

        if rve is None:
            rve = self.rve

        crd = self.nodes_coords

        all_face_xp = np.hstack(
            (
                self.face_list_xp,
                self.edge_list_xp_ym,
                self.edge_list_xp_yp,
                self.edge_list_xp_zm,
                self.edge_list_xp_zp,
                self.corner_list_xp_ym_zm,
                self.corner_list_xp_yp_zm,
                self.corner_list_xp_ym_zp,
                self.corner_list_xp_yp_zp,
            )
        )
        all_face_yp = np.hstack(
            (
                self.face_list_yp,
                self.edge_list_yp_zm,
                self.edge_list_yp_zp,
                self.edge_list_xm_yp,
                self.edge_list_xp_yp,
                self.corner_list_xm_yp_zm,
                self.corner_list_xp_yp_zm,
                self.corner_list_xm_yp_zp,
                self.corner_list_xp_yp_zp,
            )
        )
        all_face_zp = np.hstack(
            (
                self.face_list_zp,
                self.edge_list_ym_zp,
                self.edge_list_yp_zp,
                self.edge_list_xm_zp,
                self.edge_list_xp_zp,
                self.corner_list_xm_ym_zp,
                self.corner_list_xp_ym_zp,
                self.corner_list_xm_yp_zp,
                self.corner_list_xp_yp_zp,
            )
        )

        kd_trees = [
            KDTree(crd[all_face_xp]),
            KDTree(crd[all_face_yp]),
            KDTree(crd[all_face_zp]),
        ]

        offsets = [
            np.array([rve.dim_x, 0.0, 0.0]),
            np.array([0.0, rve.dim_y, 0.0]),
            np.array([0.0, 0.0, rve.dim_z]),
        ]

        faces_m = [self.face_list_xm, self.face_list_ym, self.face_list_zm]
        all_faces_p = [all_face_xp, all_face_yp, all_face_zp]

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

            dist_temp, index_temp = kd_trees[i].query(crd_face, minimum_query_points)
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
            "face_Xp": (np.asarray(index[0]), np.asarray(dist[0])),
            "face_Yp": (np.asarray(index[1]), np.asarray(dist[1])),
            "face_Zp": (np.asarray(index[2]), np.asarray(dist[2])),
        }

    def _closest_points_on_edges(
            self,
            rve: Optional[Rve] = None,
            tol: float = 1.0e-8,
    ) -> dict[str, tuple[npt.NDArray[np.int_], npt.NDArray[np.float_]]]:
        """
        Find the closest points on opposite edges to write interpolation relationship
        if a displacement condition between pair nodes is defined on such opposite surfaces
        Note : the closest neighbours is set as 2 for edges

        :param rve : RVE of the mesh bounding box. if None, the rve is built from the mesh bounding box
        """

        if rve is None:
            rve = self.rve

        crd = self.nodes_coords

        all_edge_xp_ym = np.hstack(
            (self.edge_list_xp_ym, self.corner_list_xp_ym_zm, self.corner_list_xp_ym_zp)
        )
        all_edge_xp_yp = np.hstack(
            (self.edge_list_xp_yp, self.corner_list_xp_yp_zm, self.corner_list_xp_yp_zp)
        )
        all_edge_xm_yp = np.hstack(
            (self.edge_list_xm_yp, self.corner_list_xm_yp_zm, self.corner_list_xm_yp_zp)
        )

        all_edge_xp_zm = np.hstack(
            (self.edge_list_xp_zm, self.corner_list_xp_ym_zm, self.corner_list_xp_yp_zm)
        )
        all_edge_xp_zp = np.hstack(
            (self.edge_list_xp_zp, self.corner_list_xp_ym_zp, self.corner_list_xp_yp_zp)
        )
        all_edge_xm_zp = np.hstack(
            (self.edge_list_xm_zp, self.corner_list_xm_ym_zp, self.corner_list_xm_yp_zp)
        )

        all_edge_yp_zm = np.hstack(
            (self.edge_list_ym_zm, self.corner_list_xm_yp_zm, self.corner_list_xp_yp_zm)
        )
        all_edge_yp_zp = np.hstack(
            (self.edge_list_yp_zp, self.corner_list_xm_yp_zp, self.corner_list_xp_yp_zp)
        )
        all_edge_ym_zp = np.hstack(
            (self.edge_list_ym_zp, self.corner_list_xm_ym_zp, self.corner_list_xp_ym_zp)
        )

        kd_trees = [
            KDTree(crd[all_edge_xp_ym]),
            KDTree(crd[all_edge_xp_yp]),
            KDTree(crd[all_edge_xm_yp]),
            KDTree(crd[all_edge_xp_zm]),
            KDTree(crd[all_edge_xp_zp]),
            KDTree(crd[all_edge_xm_zp]),
            KDTree(crd[all_edge_yp_zm]),
            KDTree(crd[all_edge_yp_zp]),
            KDTree(crd[all_edge_ym_zp]),
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
            self.edge_list_xm_ym,
            self.edge_list_xm_ym,
            self.edge_list_xm_ym,
            self.edge_list_xm_zm,
            self.edge_list_xm_zm,
            self.edge_list_xm_zm,
            self.edge_list_ym_zm,
            self.edge_list_ym_zm,
            self.edge_list_ym_zm,
        ]

        all_edges_p = [
            all_edge_xp_ym,
            all_edge_xp_yp,
            all_edge_xm_yp,
            all_edge_xp_zm,
            all_edge_xp_zp,
            all_edge_xm_zp,
            all_edge_yp_zm,
            all_edge_yp_zp,
            all_edge_ym_zp,
        ]

        dist = []
        index = []

        for i, edge in enumerate(edges_m):
            offset = offsets[i]
            crd_edge = crd[edge] + offset

            dist_temp, index_temp = kd_trees[i].query(crd_edge, 2)

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
            tol: float = 1.0e-8,
    ) -> dict[str, tuple[npt.NDArray[np.int_], npt.NDArray[np.float_]]]:
        """
        Find the closest points on faces and edges to write interpolation relationship
        if a displacement condition between pair nodes is defined

        :param k_neighbours : number of closest points
        :param rve : RVE of the mesh bounding box. if None, the rve is built from the mesh bounding box
        :param tol: tolerance
        """

        if isinstance(self._closest_points_on_boundaries, dict):
            return self._closest_points_on_boundaries
        if rve is None:
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
            tol: float = 1.0e-4,
    ) -> tuple[pv.PolyData, npt.NDArray[np.int_]]:
        """
        Finds boundary elements of mesh with given tolerance

        :param rve : RVE of the mesh bounding box. if None, the rve is built from the mesh bounding box
        :param tol: tolerance
        """
        if rve is None:
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
            tol: float = 1.0e-8,
    ) -> dict[str, ClosestCellsOnBoundaries]:
        """
        Find the cells to which a given point belongs to by using a ray tracing normal to the face on which it belongs
        :param rve : RVE of the mesh bounding box. if None, the rve is built from the mesh bounding box
        :param tol : tolerance to evaluate the threshold between the cells and the boundary RVE of the mesh bounding box
        """

        if rve is None:
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
            crd[self.face_list_xm],
            crd[self.face_list_ym],
            crd[self.face_list_zm],
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

            directions = np.tile(normals[i], (np.shape(face)[0], 1))

            if USE_MULTI_RAY:
                raytraceresult = surface_p.multi_ray_trace(
                    origins=face, directions=directions
                )
            else:
                intersection_points = []
                intersection_rays = []
                intersection_cells = []
                for j, face_m_i in enumerate(face):
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
