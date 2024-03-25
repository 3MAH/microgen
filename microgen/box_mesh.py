"""
Cubic mesh for FE
"""

import itertools
import warnings
from typing import Dict, List, NamedTuple, Optional, Tuple, Union

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
    closest_opposing_cells_id: List[npt.NDArray[np.int_]]


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
        elements: Dict[pv.CellType, npt.NDArray[np.int_]],
        nodes_indices: Optional[npt.NDArray[np.int_]] = None,
    ) -> None:
        super().__init__(
            nodes_coords=nodes_coords,
            elements=elements,
            nodes_indices=nodes_indices,
        )

        self.rve = self._build_rve()
        self.center: Optional[npt.NDArray[np.float_]] = None
        self._construct()

    @staticmethod
    def from_pyvista(pvmesh: Union[pv.PolyData, pv.UnstructuredGrid]):
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
        return BoxMesh(pvmesh.points, elements)

    def _construct(
        self,
        tol: float = 1.0e-8,
    ) -> None:
        """
        Construct a box Mesh with list of points in faces (excluding edges),
        edges (excluding corners) and corners.

        :param tol: tolerance to find the points on a face or an edge
        """

        crd = self.nodes_coords
        closest_point_to_rve_center = np.linalg.norm(
            crd - self.rve.center, axis=1
        ).argmin()
        self.center = crd[closest_point_to_rve_center]

        faces = {
            key + sign: np.where(np.abs(crd[:, i] - self.rve.min_point[i]) < tol)[0]
            for i, key in enumerate(["x", "y", "z"])
            for sign in ["-", "+"]
        }

        edges = {}
        for axes in itertools.combinations(["x", "y", "z"], 2):
            for signs in itertools.product(["-", "+"], repeat=2):
                ax_1 = f"{axes[0]}{signs[0]}"
                ax_2 = f"{axes[1]}{signs[1]}"
                edges[ax_1 + ax_2] = np.intersect1d(
                    faces[ax_1], faces[ax_2], assume_unique=True
                )

        self.corners = {}
        for signs in itertools.product(["-", "+"], repeat=3):
            e = f"x{signs[0]}y{signs[1]}"
            f = f"z{signs[2]}"
            self.corners[e + f] = np.intersect1d(edges[e], faces[f], assume_unique=True)

        all_corners = np.hstack(tuple(self.corners.values()))

        self.edges = {
            key: np.setdiff1d(edge, all_corners, assume_unique=True)
            for key, edge in edges.items()
        }

        all_edges_corners = np.hstack((*self.edges.values(), all_corners))

        faces = {
            key: np.setdiff1d(faces[key], all_edges_corners, assume_unique=True)
            for key in faces
        }

        for i, axes in enumerate(itertools.combinations(["x", "y", "z"], 2)):
            for signs in itertools.product(["-", "+"], repeat=2):
                key = f"{axes[0]}{signs[0]}{axes[1]}{signs[1]}"
                self.edges[key] = self.edges[key][
                    np.argsort(crd[self.edges[key], 2 - i])
                ]

        decimal_round = int(-np.log10(tol) - 1)
        self.faces = {}
        for i, axis in enumerate(["x", "y", "z"]):
            indices = [0, 1, 2]
            indices.remove(i)
            for sign in ["-", "+"]:
                key = axis + sign
                self.faces[key] = faces[key][
                    np.lexsort(
                        (
                            crd[faces[key], indices[0]],
                            crd[faces[key], indices[1]].round(decimal_round),
                        )
                    )
                ]

    def _build_rve(
        self,
    ) -> Rve:
        """
        build a representative volume element (Rve) of the considered mesh from its bounding box

        :return rve: RVE of the mesh bounding box
        """
        if not isinstance(self._pvmesh, pv.UnstructuredGrid):
            self._pvmesh = self.to_pyvista()
        xmin, xmax, ymin, ymax, zmin, zmax = self._pvmesh.bounds
        return Rve.from_min_max(
            float(xmin), float(xmax), float(ymin), float(ymax), float(zmin), float(zmax)
        )

    def _closest_points_on_faces(
        self,
        k_neighbours: int = 3,
        rve: Optional[Rve] = None,
        tol: float = 1.0e-8,
    ) -> Dict[str, Tuple[npt.NDArray[np.int_], npt.NDArray[np.float_]]]:
        """
        Find the closest points on opposite face to write interpolation relationship
        if a displacement condition between pair nodes is defined
        on such opposite surfaces
        It takes the set of points on faces 'face_Xm', 'face_Ym',
        'face_Zm' (excluding edges and corners)
        and returns for each opposite face the closest points
        (index of the points, distance) of the
        projected point on the corresponding opposite face

        :param k_neighbours : number of closest points
        :param rve : RVE of the mesh bounding box. if None,
        the rve is built from the mesh bounding box
        :param tol: tolerance

        :return dict: a dictionary with (np.array) of indices and
        a np.array of distances for each neighbor:
            'face_Xp' : (index[0], dist[0]),
            'face_Yp' : (index[1], dist[1]),
            'face_Zp' : (index[2], dist[2])
        """

        if rve is None:
            rve = self.rve

        crd = self.nodes_coords

        all_face = {
            "x+": np.hstack(
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
                )
            ),
            "y+": np.hstack(
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
                )
            ),
            "z+": np.hstack(
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
                )
            ),
        }

        kd_trees = [KDTree(crd[all_face[axis]]) for axis in ["x+", "y+", "z+"]]

        offsets = [
            np.array([rve.dim[0], 0.0, 0.0]),
            np.array([0.0, rve.dim[1], 0.0]),
            np.array([0.0, 0.0, rve.dim[2]]),
        ]

        dist = []
        index = []

        for offset, face_m, all_face_p, kd_tree in zip(
            offsets,
            [self.faces["x-"], self.faces["y-"], self.faces["z-"]],
            all_face.values(),
            kd_trees,
        ):
            crd_face = crd[face_m] + offset
            minimum_query_points = min(len(all_face_p), k_neighbours)

            if minimum_query_points < k_neighbours:
                warnings.warn(
                    "Number of query points is greater than the number of points in the KDTree."
                )

            dist_temp, index_temp = kd_tree.query(crd_face, minimum_query_points)

            index_temp_list = all_face_p[index_temp].tolist()
            if k_neighbours == 1:
                dist_temp_list = [[d] for d in dist_temp]
                index_temp_list = [[i] for i in index_temp_list]
            else:
                dist_temp_list = dist_temp.tolist()

            # If pair nodes exist (opposite node exactly match),
            # return only the pair node (the closest neighbour)
            if k_neighbours > 1:
                for i in range(len(dist_temp)):
                    if dist_temp_list[i][0] < tol:
                        dist_temp_list[i] = dist_temp_list[i][:1]
                        index_temp_list[i] = index_temp_list[i][:1]

            dist.append(dist_temp_list)
            index.append(index_temp_list)

        return {
            f"face_{axis}p": (np.asarray(index[k]), np.asarray(dist[k]))
            for k, axis in enumerate(["X", "Y", "Z"])
        }

    def _closest_points_on_edges(
        self,
        rve: Optional[Rve] = None,
        tol: float = 1.0e-8,
    ) -> Dict[str, Tuple[npt.NDArray[np.int_], npt.NDArray[np.float_]]]:
        """
        Find the closest points on opposite edges to write interpolation relationship
        if a displacement condition between pair nodes is defined on such opposite surfaces
        Note : the closest neighbours is set as 2 for edges

        :param rve : RVE of the mesh bounding box.
        if None, the rve is built from the mesh bounding box
        """

        if rve is None:
            rve = self.rve

        crd = self.nodes_coords

        all_edges = {
            "x+y-": np.hstack(
                (
                    self.edges["x+y-"],
                    self.corners["x+y-z-"],
                    self.corners["x+y-z+"],
                )
            ),
            "x+y+": np.hstack(
                (
                    self.edges["x+y+"],
                    self.corners["x+y+z-"],
                    self.corners["x+y+z+"],
                )
            ),
            "x-y+": np.hstack(
                (
                    self.edges["x-y+"],
                    self.corners["x-y+z-"],
                    self.corners["x-y+z+"],
                )
            ),
            "x+z-": np.hstack(
                (
                    self.edges["x+z-"],
                    self.corners["x+y-z-"],
                    self.corners["x+y+z-"],
                )
            ),
            "x+z+": np.hstack(
                (
                    self.edges["x+z+"],
                    self.corners["x+y-z+"],
                    self.corners["x+y+z+"],
                )
            ),
            "x-z+": np.hstack(
                (
                    self.edges["x-z+"],
                    self.corners["x-y-z+"],
                    self.corners["x-y+z+"],
                )
            ),
            "y+z-": np.hstack(
                (
                    self.edges["y+z-"],
                    self.corners["x-y+z-"],
                    self.corners["x+y+z-"],
                )
            ),
            "y+z+": np.hstack(
                (
                    self.edges["y+z+"],
                    self.corners["x-y+z+"],
                    self.corners["x+y+z+"],
                )
            ),
            "y-z+": np.hstack(
                (
                    self.edges["y-z+"],
                    self.corners["x-y-z+"],
                    self.corners["x+y-z+"],
                )
            ),
        }

        kd_trees = [KDTree(crd[edge]) for edge in all_edges.values()]

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
            self.edges[f"{ax[0]}-{ax[1]}-"]
            for ax in itertools.combinations(("x", "y", "z"), 2)
            for _ in range(3)
        ]

        all_edges_p = [
            all_edges[f"{ax[0]}{d[0]}{ax[1]}{d[1]}"]
            for ax in itertools.combinations(("x", "y", "z"), 2)
            for d in itertools.product(("-", "+"), repeat=2)
            if d != ("-", "-")
        ]

        dist = []
        index = []
        for edge, offset, kd_tree, all_edge_p in zip(
            edges_m, offsets, kd_trees, all_edges_p
        ):
            dist_temp, index_temp = kd_tree.query(crd[edge] + offset, 2)

            dist_temp_list = dist_temp.tolist()
            index_temp_list = all_edge_p[index_temp].tolist()
            for j, (dist_tmp, index_tmp) in enumerate(
                zip(dist_temp_list, index_temp_list)
            ):
                if dist_temp_list[j][0] < tol:
                    dist_temp_list[j] = dist_tmp[:1]
                    index_temp_list[j] = index_tmp[:1]

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
    ) -> Dict[str, Tuple[npt.NDArray[np.int_], npt.NDArray[np.float_]]]:
        """
        Find the closest points on faces and edges to write interpolation relationship
        if a displacement condition between pair nodes is defined

        :param k_neighbours : number of closest points
        :param rve : RVE of the mesh bounding box.
        if None, the rve is built from the mesh bounding box
        :param tol: tolerance
        """

        if rve is None:
            rve = self.rve

        dict_faces = self._closest_points_on_faces(k_neighbours, rve, tol)
        dict_edges = self._closest_points_on_edges(rve, tol)

        return {**dict_faces, **dict_edges}

    def boundary_elements(
        self,
        rve: Optional[Rve] = None,
        tol: float = 1.0e-4,
    ) -> Tuple[pv.PolyData, npt.NDArray[np.int_]]:
        """
        Finds boundary elements of mesh with given tolerance

        :param rve : RVE of the mesh bounding box.
        if None, the rve is built from the mesh bounding box
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
            np.array([rve.min_point[0], rve.center[1], rve.center[2]]),
            np.array([rve.max_point[0], rve.center[1], rve.center[2]]),
            np.array([rve.center[0], rve.min_point[1], rve.center[2]]),
            np.array([rve.center[0], rve.max_point[1], rve.center[2]]),
            np.array([rve.center[0], rve.center[1], rve.min_point[2]]),
            np.array([rve.center[0], rve.center[1], rve.max_point[2]]),
        ]
        size_planes = [
            2.0 * rve.dim[0],
            2.0 * rve.dim[0],
            2.0 * rve.dim[1],
            2.0 * rve.dim[1],
            2.0 * rve.dim[2],
            2.0 * rve.dim[2],
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
    ) -> Dict[str, ClosestCellsOnBoundaries]:
        """
        Find the cells to which a given point belongs to by using a ray tracing
        normal to the face on which it belongs
        :param rve : RVE of the mesh bounding box. if None, the rve is built
        from the mesh bounding box
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
            np.array([rve.max_point[0], rve.center[1], rve.center[2]]),
            np.array([rve.center[0], rve.max_point[1], rve.center[2]]),
            np.array([rve.center[0], rve.center[1], rve.max_point[2]]),
        ]
        size_planes = [2.0 * rve.dim[0], 2.0 * rve.dim[1], 2.0 * rve.dim[2]]

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
            f"face_{axis}m": ClosestCellsOnBoundaries(
                list_cells_for_each_face[i],
                list_ray_trace[i][0],
                list_cells_for_each_face[i][list_ray_trace[i][2]],
            )
            for i, axis in enumerate(["X", "Y", "Z"])
        }
