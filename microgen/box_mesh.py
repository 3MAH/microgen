"""Cubic mesh for FE."""

from __future__ import annotations

import warnings
from itertools import product
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

AXES = ("x", "y", "z")


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

    def _construct(self: BoxMesh, tol: float = 1.0e-8) -> None:
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

        faces = {
            f"{axis}{sign}": np.where(
                np.abs(crd[:, a] - bound[a]) < tol,
            )[0]
            for a, axis in enumerate(("x", "y", "z"))
            for bound, sign in zip((self.rve.min_point, self.rve.max_point), ("-", "+"))
        }

        edges = {
            f"{axis1}{sign1}{axis2}{sign2}": np.intersect1d(
                faces[f"{axis1}{sign1}"],
                faces[f"{axis2}{sign2}"],
                assume_unique=True,
            )
            for (axis1, axis2) in [("x", "y"), ("x", "z"), ("y", "z")]  # axis pairs
            for (sign1, sign2) in product(("-", "+"), repeat=2)  # sign combinations
        }

        # extract corners from the intersection of edges
        corners = {
            f"x{sx}y{sy}z{sz}": np.intersect1d(
                edges[f"x{sx}y{sy}"],
                edges[f"y{sy}z{sz}"],
                assume_unique=True,
            )
            for sx, sy, sz in product(("-", "+"), repeat=3)
        }

        # Remove nodes that belong to several sets
        all_corners = np.hstack(list(corners.values()))

        for axis1, axis2 in [("x", "y"), ("x", "z"), ("y", "z")]:
            for sign1, sign2 in product(("-", "+"), repeat=2):
                edges[f"{axis1}{sign1}{axis2}{sign2}"] = np.setdiff1d(
                    edges[f"{axis1}{sign1}{axis2}{sign2}"],
                    all_corners,
                    assume_unique=True,
                )

        all_edges_corners = np.hstack([*list(edges.values()), all_corners])

        decimal_round = int(-np.log10(tol) - 1)
        for a, axis in enumerate(["x", "y", "z"]):
            for sign in ["-", "+"]:
                faces[f"{axis}{sign}"] = np.setdiff1d(
                    faces[f"{axis}{sign}"],
                    all_edges_corners,
                    assume_unique=True,
                )
                faces[f"{axis}{sign}"] = faces[f"{axis}{sign}"][
                    np.lexsort(
                        (
                            crd[faces[f"{axis}{sign}"], (a + 1) % 3],
                            crd[faces[f"{axis}{sign}"], (a + 2) % 3].round(
                                decimal_round,
                            ),
                        ),
                    )
                ]

        self.corners = corners
        self.edges = edges
        self.faces = faces

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
            'x+': (index[0], dist[0]),
            'y+': (index[1], dist[1]),
            'z+': (index[2], dist[2])
        """
        if rve is None:
            rve = self.rve

        crd = self.nodes_coords

        all_face_xp = np.hstack(
            (
                self.faces["x+"],
                self.edges["x+y-"],
                self.edges["x+y+"],
                self.edges["x+z-"],
                self.edges["x+z+"],
                self.corners["x+y-z-"],
                self.corners["x+y+z-"],
                self.corners["x+y-z+"],
                self.corners["x+y+z+"],
            ),
        )
        all_face_yp = np.hstack(
            (
                self.faces["y+"],
                self.edges["y+z-"],
                self.edges["y+z+"],
                self.edges["x-y+"],
                self.edges["x+y+"],
                self.corners["x-y+z-"],
                self.corners["x+y+z-"],
                self.corners["x-y+z+"],
                self.corners["x+y+z+"],
            ),
        )
        all_face_zp = np.hstack(
            (
                self.faces["z+"],
                self.edges["y-z+"],
                self.edges["y+z+"],
                self.edges["x-z+"],
                self.edges["x+z+"],
                self.corners["x-y-z+"],
                self.corners["x+y-z+"],
                self.corners["x-y+z+"],
                self.corners["x+y+z+"],
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

        faces_m = [self.faces["x-"], self.faces["y-"], self.faces["z-"]]
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
            "x+": (np.asarray(index[0]), np.asarray(dist[0])),
            "y+": (np.asarray(index[1]), np.asarray(dist[1])),
            "z+": (np.asarray(index[2]), np.asarray(dist[2])),
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
                self.edges["x+y-"],
                self.corners["x+y-z-"],
                self.corners["x+y-z+"],
            ),
        )
        all_edge_xp_yp = np.hstack(
            (
                self.edges["x+y+"],
                self.corners["x+y+z-"],
                self.corners["x+y+z+"],
            ),
        )
        all_edge_xm_yp = np.hstack(
            (
                self.edges["x-y+"],
                self.corners["x-y+z-"],
                self.corners["x-y+z+"],
            ),
        )

        all_edge_xp_zm = np.hstack(
            (
                self.edges["x+z-"],
                self.corners["x+y-z-"],
                self.corners["x+y+z-"],
            ),
        )
        all_edge_xp_zp = np.hstack(
            (
                self.edges["x+z+"],
                self.corners["x+y-z+"],
                self.corners["x+y+z+"],
            ),
        )
        all_edge_xm_zp = np.hstack(
            (
                self.edges["x-z+"],
                self.corners["x-y-z+"],
                self.corners["x-y+z+"],
            ),
        )

        all_edge_yp_zm = np.hstack(
            (
                self.edges["y+z-"],
                self.corners["x-y+z-"],
                self.corners["x+y+z-"],
            ),
        )
        all_edge_yp_zp = np.hstack(
            (
                self.edges["y+z+"],
                self.corners["x-y+z+"],
                self.corners["x+y+z+"],
            ),
        )
        all_edge_ym_zp = np.hstack(
            (
                self.edges["y-z+"],
                self.corners["x-y-z+"],
                self.corners["x+y-z+"],
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
            self.edges["x-y-"],
            self.edges["x-y-"],
            self.edges["x-y-"],
            self.edges["x-z-"],
            self.edges["x-z-"],
            self.edges["x-z-"],
            self.edges["y-z-"],
            self.edges["y-z-"],
            self.edges["y-z-"],
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
            "x+y-": (index[0], dist[0]),
            "x+y+": (index[1], dist[1]),
            "x-y+": (index[2], dist[2]),
            "x+z-": (index[3], dist[3]),
            "x+z+": (index[4], dist[4]),
            "x-z+": (index[5], dist[5]),
            "y+z-": (index[6], dist[6]),
            "y+z+": (index[7], dist[7]),
            "y-z+": (index[8], dist[8]),
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

        normals = np.eye(3)
        origins_p = np.array(
            [
                [rve.max_point[i] if i == j else rve.center[i] for j in range(3)]
                for i in range(3)
            ],
        )
        size_planes = 2.0 * rve.dim

        faces_m = [crd[self.faces[f"{axis}-"]] for axis in AXES]
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
            f"{axis}+": ClosestCellsOnBoundaries(
                list_cells_for_each_face[a],
                list_ray_trace[a][0],
                list_cells_for_each_face[a][list_ray_trace[a][2]],
            )
            for a, axis in enumerate(AXES)
        }
