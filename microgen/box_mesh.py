"""Cubic mesh for FE."""

from __future__ import annotations

import warnings
from typing import NamedTuple

import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.spatial import KDTree

from .rve import Rve
from .single_mesh import SingleMesh, check_if_only_linear_tetrahedral

# ruff: noqa: ERA001
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
    """Class to manage closest cells on boundaries.

    :param cells_id: list of cells identification number
    :param intersection_point_coords: list of intersection points of projected
    rays from one face's nodes with the opposing face
    :param closest_opposing_cells_id:
        list of closest cells identification number on opposing face
    """

    cells_id: npt.NDArray[np.int_]
    intersection_point_coords: npt.NDArray[np.float64]
    closest_opposing_cells_id: list[npt.NDArray[np.int_]]


class BoxMesh(SingleMesh):
    """Class to manage list of Nodes and Elements inside a Rve.

    :param nodes_coords: list of nodes (np.ndarray)
    :param elements: list of elements (np.ndarray)
    :param nodes_indices : index of node list
        (if different from the natural index of nodes array)
    """

    def __init__(
        self: BoxMesh,
        nodes_coords: npt.NDArray[np.float64],
        elements: dict[pv.CellType, npt.NDArray[np.int_]],
        nodes_indices: npt.NDArray[np.int_] | None = None,
    ) -> None:
        """Initialize the BoxMesh object with nodes and elements."""
        super().__init__(
            nodes_coords=nodes_coords,
            elements=elements,
            nodes_indices=nodes_indices,
        )

        self.rve = self._build_rve()
        self.center: npt.NDArray[np.float64] | None = None
        self._construct()

    @staticmethod
    def from_pyvista(pvmesh: pv.PolyData | pv.UnstructuredGrid) -> BoxMesh:
        """Build a BoxMesh from a pyvista UnstructuredGrid or PolyData mesh.

        If pvmesh is a PolyData object, it is cast into an Unstructured Grid).
        Node and element data are not copied (by default a shallow copy is operated)
        Mesh with multiple element type are not handled.
        Note that for now singleMesh works only considering tetrahedral elements

        :param pvmesh: the mesh as a pyvista UnstructuredGrid object

        :return : A BoxMesh Object from the tetrahedral elements
            of the pyvista Unstructured Grid
        """
        if isinstance(pvmesh, pv.PolyData):
            pvmesh = pvmesh.cast_to_unstructured_grid()

        # extract only the tetrahedral elements
        # raises an Exception if 3d elements other than
        # linear tetrahedra are found in the mesh
        check_if_only_linear_tetrahedral(pvmesh)
        elements = {pv.CellType.TETRA: pvmesh.cells_dict[pv.CellType.TETRA]}
        return BoxMesh(pvmesh.points, elements)

    def _construct(self: BoxMesh, tol: float = 1.0e-8) -> None:  # noqa: PLR0915 TODO(kmarchais): refactor
        """Construct a box Mesh.

        With list of points in faces (excluding edges),
        edges (excluding corners) and corners.

        :param tol: tolerance to find the points on a face or an edge
        """
        crd = self.nodes_coords
        closest_point_to_rve_center = np.linalg.norm(
            crd - self.rve.center,
            axis=1,
        ).argmin()
        self.center = crd[closest_point_to_rve_center]

        face_list_xm = np.where(np.abs(crd[:, 0] - self.rve.min_point[0]) < tol)[0]
        face_list_xp = np.where(np.abs(crd[:, 0] - self.rve.max_point[0]) < tol)[0]

        face_list_ym = np.where(np.abs(crd[:, 1] - self.rve.min_point[1]) < tol)[0]
        face_list_yp = np.where(np.abs(crd[:, 1] - self.rve.max_point[1]) < tol)[0]

        face_list_zm = np.where(np.abs(crd[:, 2] - self.rve.min_point[2]) < tol)[0]
        face_list_zp = np.where(np.abs(crd[:, 2] - self.rve.max_point[2]) < tol)[0]

        edge_list_xm_ym = np.intersect1d(face_list_xm, face_list_ym, assume_unique=True)
        edge_list_xp_ym = np.intersect1d(face_list_xp, face_list_ym, assume_unique=True)
        edge_list_xp_yp = np.intersect1d(face_list_xp, face_list_yp, assume_unique=True)
        edge_list_xm_yp = np.intersect1d(face_list_xm, face_list_yp, assume_unique=True)

        edge_list_xm_zm = np.intersect1d(face_list_xm, face_list_zm, assume_unique=True)
        edge_list_xp_zm = np.intersect1d(face_list_xp, face_list_zm, assume_unique=True)
        edge_list_xp_zp = np.intersect1d(face_list_xp, face_list_zp, assume_unique=True)
        edge_list_xm_zp = np.intersect1d(face_list_xm, face_list_zp, assume_unique=True)

        edge_list_ym_zm = np.intersect1d(face_list_ym, face_list_zm, assume_unique=True)
        edge_list_yp_zm = np.intersect1d(face_list_yp, face_list_zm, assume_unique=True)
        edge_list_yp_zp = np.intersect1d(face_list_yp, face_list_zp, assume_unique=True)
        edge_list_ym_zp = np.intersect1d(face_list_ym, face_list_zp, assume_unique=True)

        # extract corners from the intersection of edges
        corner_list_xm_ym_zm = np.intersect1d(
            edge_list_xm_ym,
            edge_list_ym_zm,
            assume_unique=True,
        )
        corner_list_xm_yp_zm = np.intersect1d(
            edge_list_xm_yp,
            edge_list_yp_zm,
            assume_unique=True,
        )
        corner_list_xp_ym_zm = np.intersect1d(
            edge_list_xp_ym,
            edge_list_ym_zm,
            assume_unique=True,
        )
        corner_list_xp_yp_zm = np.intersect1d(
            edge_list_xp_yp,
            edge_list_yp_zm,
            assume_unique=True,
        )
        corner_list_xm_ym_zp = np.intersect1d(
            edge_list_xm_ym,
            edge_list_ym_zp,
            assume_unique=True,
        )
        corner_list_xm_yp_zp = np.intersect1d(
            edge_list_xm_yp,
            edge_list_yp_zp,
            assume_unique=True,
        )
        corner_list_xp_ym_zp = np.intersect1d(
            edge_list_xp_ym,
            edge_list_ym_zp,
            assume_unique=True,
        )
        corner_list_xp_yp_zp = np.intersect1d(
            edge_list_xp_yp,
            edge_list_yp_zp,
            assume_unique=True,
        )

        # Remove nodes that belong to several sets
        all_corners = np.hstack(
            (
                corner_list_xm_ym_zm,
                corner_list_xm_yp_zm,
                corner_list_xp_ym_zm,
                corner_list_xp_yp_zm,
                corner_list_xm_ym_zp,
                corner_list_xm_yp_zp,
                corner_list_xp_ym_zp,
                corner_list_xp_yp_zp,
            ),
        )

        edge_list_xm_ym = np.setdiff1d(edge_list_xm_ym, all_corners, assume_unique=True)
        edge_list_xp_ym = np.setdiff1d(edge_list_xp_ym, all_corners, assume_unique=True)
        edge_list_xp_yp = np.setdiff1d(edge_list_xp_yp, all_corners, assume_unique=True)
        edge_list_xm_yp = np.setdiff1d(edge_list_xm_yp, all_corners, assume_unique=True)

        edge_list_xm_zm = np.setdiff1d(edge_list_xm_zm, all_corners, assume_unique=True)
        edge_list_xp_zm = np.setdiff1d(edge_list_xp_zm, all_corners, assume_unique=True)
        edge_list_xp_zp = np.setdiff1d(edge_list_xp_zp, all_corners, assume_unique=True)
        edge_list_xm_zp = np.setdiff1d(edge_list_xm_zp, all_corners, assume_unique=True)

        edge_list_ym_zm = np.setdiff1d(edge_list_ym_zm, all_corners, assume_unique=True)
        edge_list_yp_zm = np.setdiff1d(edge_list_yp_zm, all_corners, assume_unique=True)
        edge_list_yp_zp = np.setdiff1d(edge_list_yp_zp, all_corners, assume_unique=True)
        edge_list_ym_zp = np.setdiff1d(edge_list_ym_zp, all_corners, assume_unique=True)

        all_edges_corners = np.hstack(
            (
                edge_list_xm_ym,
                edge_list_xp_ym,
                edge_list_xp_yp,
                edge_list_xm_yp,
                edge_list_xm_zm,
                edge_list_xp_zm,
                edge_list_xp_zp,
                edge_list_xm_zp,
                edge_list_ym_zm,
                edge_list_yp_zm,
                edge_list_yp_zp,
                edge_list_ym_zp,
                all_corners,
            ),
        )

        face_list_xm = np.setdiff1d(face_list_xm, all_edges_corners, assume_unique=True)
        face_list_xp = np.setdiff1d(face_list_xp, all_edges_corners, assume_unique=True)
        face_list_ym = np.setdiff1d(face_list_ym, all_edges_corners, assume_unique=True)
        face_list_yp = np.setdiff1d(face_list_yp, all_edges_corners, assume_unique=True)
        face_list_zm = np.setdiff1d(face_list_zm, all_edges_corners, assume_unique=True)
        face_list_zp = np.setdiff1d(face_list_zp, all_edges_corners, assume_unique=True)

        edge_list_xm_ym = edge_list_xm_ym[np.argsort(crd[edge_list_xm_ym, 2])]
        edge_list_xp_ym = edge_list_xp_ym[np.argsort(crd[edge_list_xp_ym, 2])]
        edge_list_xp_yp = edge_list_xp_yp[np.argsort(crd[edge_list_xp_yp, 2])]
        edge_list_xm_yp = edge_list_xm_yp[np.argsort(crd[edge_list_xm_yp, 2])]

        edge_list_xm_zm = edge_list_xm_zm[np.argsort(crd[edge_list_xm_zm, 1])]
        edge_list_xp_zm = edge_list_xp_zm[np.argsort(crd[edge_list_xp_zm, 1])]
        edge_list_xp_zp = edge_list_xp_zp[np.argsort(crd[edge_list_xp_zp, 1])]
        edge_list_xm_zp = edge_list_xm_zp[np.argsort(crd[edge_list_xm_zp, 1])]

        edge_list_ym_zm = edge_list_ym_zm[np.argsort(crd[edge_list_ym_zm, 0])]
        edge_list_yp_zm = edge_list_yp_zm[np.argsort(crd[edge_list_yp_zm, 0])]
        edge_list_yp_zp = edge_list_yp_zp[np.argsort(crd[edge_list_yp_zp, 0])]
        edge_list_ym_zp = edge_list_ym_zp[np.argsort(crd[edge_list_ym_zp, 0])]

        decimal_round = int(-np.log10(tol) - 1)
        face_list_xm = face_list_xm[
            np.lexsort(
                (
                    crd[face_list_xm, 1],
                    crd[face_list_xm, 2].round(decimal_round),
                ),
            )
        ]
        face_list_xp = face_list_xp[
            np.lexsort(
                (
                    crd[face_list_xp, 1],
                    crd[face_list_xp, 2].round(decimal_round),
                ),
            )
        ]
        face_list_ym = face_list_ym[
            np.lexsort(
                (
                    crd[face_list_ym, 0],
                    crd[face_list_ym, 2].round(decimal_round),
                ),
            )
        ]
        face_list_yp = face_list_yp[
            np.lexsort(
                (
                    crd[face_list_yp, 0],
                    crd[face_list_yp, 2].round(decimal_round),
                ),
            )
        ]
        face_list_zm = face_list_zm[
            np.lexsort(
                (
                    crd[face_list_zm, 0],
                    crd[face_list_zm, 1].round(decimal_round),
                ),
            )
        ]
        face_list_zp = face_list_zp[
            np.lexsort(
                (
                    crd[face_list_zp, 0],
                    crd[face_list_zp, 1].round(decimal_round),
                ),
            )
        ]

        self.corners = {
            "corner_xm_ym_zm": corner_list_xm_ym_zm,
            "corner_xm_yp_zm": corner_list_xm_yp_zm,
            "corner_xp_ym_zm": corner_list_xp_ym_zm,
            "corner_xp_yp_zm": corner_list_xp_yp_zm,
            "corner_xm_ym_zp": corner_list_xm_ym_zp,
            "corner_xm_yp_zp": corner_list_xm_yp_zp,
            "corner_xp_ym_zp": corner_list_xp_ym_zp,
            "corner_xp_yp_zp": corner_list_xp_yp_zp,
        }

        self.edges = {
            "edge_xm_ym": edge_list_xm_ym,
            "edge_xp_ym": edge_list_xp_ym,
            "edge_xp_yp": edge_list_xp_yp,
            "edge_xm_yp": edge_list_xm_yp,
            "edge_xm_zm": edge_list_xm_zm,
            "edge_xp_zm": edge_list_xp_zm,
            "edge_xp_zp": edge_list_xp_zp,
            "edge_xm_zp": edge_list_xm_zp,
            "edge_ym_zm": edge_list_ym_zm,
            "edge_yp_zm": edge_list_yp_zm,
            "edge_yp_zp": edge_list_yp_zp,
            "edge_ym_zp": edge_list_ym_zp,
        }

        self.faces = {
            "face_xm": face_list_xm,
            "face_xp": face_list_xp,
            "face_ym": face_list_ym,
            "face_yp": face_list_yp,
            "face_zm": face_list_zm,
            "face_zp": face_list_zp,
        }

    def _build_rve(self: BoxMesh) -> Rve:
        """Build a representative volume element (Rve) from the mesh's bounding box.

        :return rve: RVE of the mesh bounding box
        """
        if not isinstance(self._pvmesh, pv.UnstructuredGrid):
            self._pvmesh = self.to_pyvista()
        xmin, xmax, ymin, ymax, zmin, zmax = self._pvmesh.bounds
        return Rve.from_min_max(
            float(xmin),
            float(xmax),
            float(ymin),
            float(ymax),
            float(zmin),
            float(zmax),
        )

    def _closest_points_on_faces(
        self: BoxMesh,
        k_neighbours: int = 3,
        rve: Rve | None = None,
        tol: float = 1.0e-8,
    ) -> dict[str, tuple[npt.NDArray[np.int_], npt.NDArray[np.float64]]]:
        """Find closest points on opposite faces to write interpolation relationship.

        If a displacement condition between pair nodes is defined on such opposite
        surfaces, it takes the set of points on faces 'face_Xm', 'face_Ym', 'face_Zm'
        (excluding edges and corners) and returns for each opposite face the closest
        points (index of the points, distance) of the projected point on the
        corresponding opposite face.

        :param k_neighbours: Number of closest points.
        :param rve: RVE of the mesh bounding box.
            If None, the RVE is built from the mesh bounding box.
        :param tol: Tolerance.

        :return dict:A dictionary with (np.array) of indices and
            a np.array of distances for each neighbor:
            'face_xp': (index[0], dist[0]),
            'face_xp': (index[1], dist[1]),
            'face_xp': (index[2], dist[2])
        """
        if rve is None:
            rve = self.rve

        crd = self.nodes_coords

        all_face_xp = np.hstack(
            (
                self.faces["face_xp"],
                self.edges["edge_xp_ym"],
                self.edges["edge_xp_yp"],
                self.edges["edge_xp_zm"],
                self.edges["edge_xp_zp"],
                self.corners["corner_xp_ym_zm"],
                self.corners["corner_xp_yp_zm"],
                self.corners["corner_xp_ym_zp"],
                self.corners["corner_xp_yp_zp"],
            ),
        )
        all_face_yp = np.hstack(
            (
                self.faces["face_yp"],
                self.edges["edge_yp_zm"],
                self.edges["edge_yp_zp"],
                self.edges["edge_xm_yp"],
                self.edges["edge_xp_yp"],
                self.corners["corner_xm_yp_zm"],
                self.corners["corner_xp_yp_zm"],
                self.corners["corner_xm_yp_zp"],
                self.corners["corner_xp_yp_zp"],
            ),
        )
        all_face_zp = np.hstack(
            (
                self.faces["face_zp"],
                self.edges["edge_ym_zp"],
                self.edges["edge_yp_zp"],
                self.edges["edge_xm_zp"],
                self.edges["edge_xp_zp"],
                self.corners["corner_xm_ym_zp"],
                self.corners["corner_xp_ym_zp"],
                self.corners["corner_xm_yp_zp"],
                self.corners["corner_xp_yp_zp"],
            ),
        )

        kd_trees = [
            KDTree(crd[all_face_xp]),
            KDTree(crd[all_face_yp]),
            KDTree(crd[all_face_zp]),
        ]

        offsets = [
            np.array([rve.dim[0], 0.0, 0.0]),
            np.array([0.0, rve.dim[1], 0.0]),
            np.array([0.0, 0.0, rve.dim[2]]),
        ]

        faces_m = [self.faces["face_xm"], self.faces["face_ym"], self.faces["face_zm"]]
        all_faces_p = [all_face_xp, all_face_yp, all_face_zp]

        dist = []
        index = []

        for i, face in enumerate(faces_m):
            offset = offsets[i]
            crd_face = crd[face] + offset
            minimum_query_points = min(len(all_faces_p[i]), k_neighbours)

            if minimum_query_points < k_neighbours:
                warnings.warn(
                    (
                        "Number of query points is greater than the number "
                        "of points in the KDTree."
                    ),
                    stacklevel=2,
                )

            dist_temp, index_temp = kd_trees[i].query(crd_face, minimum_query_points)
            all_faces_p_i = all_faces_p[i]
            if k_neighbours == 1:
                dist_temp_list = list(dist_temp)
                index_temp_list = [all_faces_p_i[index_temp]]
            else:
                dist_temp_list = dist_temp.tolist()
                index_temp_list = all_faces_p_i[index_temp].tolist()

            # If pair nodes exist (opposite node exactly match)
            # return only the pair node (the closest neighbour)
            if k_neighbours > 1:
                for nb_pts_on_face in range(len(dist_temp)):
                    if dist_temp_list[nb_pts_on_face][0] < tol:
                        dist_temp_list[nb_pts_on_face] = dist_temp_list[nb_pts_on_face][
                            :1
                        ]
                        index_temp_list[nb_pts_on_face] = index_temp_list[
                            nb_pts_on_face
                        ][:1]

            dist.append(dist_temp_list)
            index.append(index_temp_list)

        return {
            "face_xp": (np.asarray(index[0]), np.asarray(dist[0])),
            "face_yp": (np.asarray(index[1]), np.asarray(dist[1])),
            "face_zp": (np.asarray(index[2]), np.asarray(dist[2])),
        }

    def _closest_points_on_edges(
        self: BoxMesh,
        rve: Rve | None = None,
        tol: float = 1.0e-8,
    ) -> dict[str, tuple[npt.NDArray[np.int_], npt.NDArray[np.float64]]]:
        """Find closest points on opposite edges to write interpolation relationship.

        if a displacement condition between pair nodes is
        defined on such opposite surfaces
        Note : the closest neighbours is set as 2 for edges

        :param rve : RVE of the mesh bounding box.
            if None, the rve is built from the mesh bounding box
        """
        if rve is None:
            rve = self.rve

        crd = self.nodes_coords

        all_edge_xp_ym = np.hstack(
            (
                self.edges["edge_xp_ym"],
                self.corners["corner_xp_ym_zm"],
                self.corners["corner_xp_ym_zp"],
            ),
        )
        all_edge_xp_yp = np.hstack(
            (
                self.edges["edge_xp_yp"],
                self.corners["corner_xp_yp_zm"],
                self.corners["corner_xp_yp_zp"],
            ),
        )
        all_edge_xm_yp = np.hstack(
            (
                self.edges["edge_xm_yp"],
                self.corners["corner_xm_yp_zm"],
                self.corners["corner_xm_yp_zp"],
            ),
        )

        all_edge_xp_zm = np.hstack(
            (
                self.edges["edge_xp_zm"],
                self.corners["corner_xp_ym_zm"],
                self.corners["corner_xp_yp_zm"],
            ),
        )
        all_edge_xp_zp = np.hstack(
            (
                self.edges["edge_xp_zp"],
                self.corners["corner_xp_ym_zp"],
                self.corners["corner_xp_yp_zp"],
            ),
        )
        all_edge_xm_zp = np.hstack(
            (
                self.edges["edge_xm_zp"],
                self.corners["corner_xm_ym_zp"],
                self.corners["corner_xm_yp_zp"],
            ),
        )

        all_edge_yp_zm = np.hstack(
            (
                self.edges["edge_yp_zm"],
                self.corners["corner_xm_yp_zm"],
                self.corners["corner_xp_yp_zm"],
            ),
        )
        all_edge_yp_zp = np.hstack(
            (
                self.edges["edge_yp_zp"],
                self.corners["corner_xm_yp_zp"],
                self.corners["corner_xp_yp_zp"],
            ),
        )
        all_edge_ym_zp = np.hstack(
            (
                self.edges["edge_ym_zp"],
                self.corners["corner_xm_ym_zp"],
                self.corners["corner_xp_ym_zp"],
            ),
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
            np.array([rve.dim[0], 0.0, 0.0]),
            np.array([rve.dim[0], rve.dim[1], 0.0]),
            np.array([0.0, rve.dim[1], 0.0]),
            np.array([rve.dim[0], 0.0, 0.0]),
            np.array([rve.dim[0], 0.0, rve.dim[2]]),
            np.array([0.0, 0.0, rve.dim[2]]),
            np.array([0.0, rve.dim[1], 0.0]),
            np.array([0.0, rve.dim[1], rve.dim[2]]),
            np.array([0.0, 0.0, rve.dim[2]]),
        ]

        edges_m = [
            self.edges["edge_xm_ym"],
            self.edges["edge_xm_ym"],
            self.edges["edge_xm_ym"],
            self.edges["edge_xm_zm"],
            self.edges["edge_xm_zm"],
            self.edges["edge_xm_zm"],
            self.edges["edge_ym_zm"],
            self.edges["edge_ym_zm"],
            self.edges["edge_ym_zm"],
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

            dist_tmp, index_tmp = kd_trees[i].query(crd_edge, 2)

            dist_tmp_list = dist_tmp.tolist()
            index_tmp_list = all_edges_p[i][index_tmp].tolist()

            for nb_pts_on_edge in range(len(dist_tmp)):
                if dist_tmp_list[nb_pts_on_edge][0] < tol:
                    dist_tmp_list[nb_pts_on_edge] = dist_tmp_list[nb_pts_on_edge][:1]
                    index_tmp_list[nb_pts_on_edge] = index_tmp_list[nb_pts_on_edge][:1]

            dist.append(dist_tmp_list)
            index.append(index_tmp_list)

        return {
            "edge_xpym": (index[0], dist[0]),
            "edge_xpyp": (index[1], dist[1]),
            "edge_xmyp": (index[2], dist[2]),
            "edge_xpzm": (index[3], dist[3]),
            "edge_xpzp": (index[4], dist[4]),
            "edge_xmzp": (index[5], dist[5]),
            "edge_ypzm": (index[6], dist[6]),
            "edge_ypzp": (index[7], dist[7]),
            "edge_ymzp": (index[8], dist[8]),
        }

    def closest_points_on_boundaries(
        self: BoxMesh,
        k_neighbours: int = 3,
        rve: Rve | None = None,
        tol: float = 1.0e-8,
    ) -> dict[str, tuple[npt.NDArray[np.int_], npt.NDArray[np.float64]]]:
        """Find closest points on faces and edges to write interpolation relationship.

        if a displacement condition between pair nodes is defined

        :param k_neighbours : number of closest points
        :param rve : RVE of the mesh bounding box.
            If None, the rve is built from the mesh bounding box
        :param tol: tolerance
        """
        if rve is None:
            rve = self.rve

        dict_faces = self._closest_points_on_faces(k_neighbours, rve, tol)
        dict_edges = self._closest_points_on_edges(rve, tol)

        return {**dict_faces, **dict_edges}

    def boundary_elements(
        self: BoxMesh,
        rve: Rve | None = None,
        tol: float = 1.0e-4,
    ) -> tuple[pv.PolyData, npt.NDArray[np.int_]]:
        """Find boundary elements of mesh with given tolerance.

        :param rve : RVE of the mesh bounding box.
            If None, the rve is built from the mesh bounding box
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
            2.0 * np.max([rve.dim[1], rve.dim[2]]),
            2.0 * np.max([rve.dim[1], rve.dim[2]]),
            2.0 * np.max([rve.dim[0], rve.dim[2]]),
            2.0 * np.max([rve.dim[0], rve.dim[2]]),
            2.0 * np.max([rve.dim[0], rve.dim[1]]),
            2.0 * np.max([rve.dim[0], rve.dim[1]]),
        ]

        surface = self.surface
        surface["CellIDs"] = np.arange(surface.n_cells)

        boundary_elements = pv.PolyData()

        for normal, origin_p, size_plane in zip(normals, origins_p, size_planes):
            plane = pv.Plane(
                center=origin_p,
                direction=normal,
                i_size=size_plane,
                j_size=size_plane,
            )
            surface.compute_implicit_distance(plane, inplace=True)
            surface_p = surface.threshold(
                [-tol, tol],
                all_scalars=True,
                scalars="implicit_distance",
            ).extract_surface()
            boundary_elements = boundary_elements.append_polydata(surface_p)

        return boundary_elements, boundary_elements["CellIDs"]

    def closest_cells_on_boundaries(
        self: BoxMesh,
        rve: Rve | None = None,
        tol: float = 1.0e-8,
    ) -> dict[str, ClosestCellsOnBoundaries]:
        """Find the cells containing a point using a ray tracing normal to its face.

        :param rve : RVE of the mesh bounding box.
            If None, the rve is built from the mesh bounding box
        :param tol : tolerance to evaluate the threshold between the cells and
            the boundary RVE of the mesh bounding box
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
            crd[self.faces["face_xm"]],
            crd[self.faces["face_ym"]],
            crd[self.faces["face_zm"]],
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
                [-tol, tol],
                all_scalars=True,
                scalars="implicit_distance",
            )
            list_cells_for_each_face.append(surface_p["CellIDs"])

            directions = np.tile(normals[i], (np.shape(face)[0], 1))

            if USE_MULTI_RAY:
                raytraceresult = surface_p.multi_ray_trace(
                    origins=face,
                    directions=directions,
                )
            else:
                intersection_points = []
                intersection_rays = []
                intersection_cells = []
                for j, face_m_i in enumerate(face):
                    end_point = face_m_i + size_planes[i] * directions[j]
                    intersection_point, intersection_cell = surface_p.ray_trace(
                        origin=face_m_i,
                        end_point=end_point,
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
