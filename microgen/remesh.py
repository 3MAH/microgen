import gmsh
from microgen import Rve
import numpy as np
import math as m
from collections import namedtuple

Triangle = namedtuple("Triangle", ["node1", "node2", "node3", "tag"])

_MESH_DIM = 3
_BOUNDARY_DIM = 2


def remesh_periodic(
    mesh_file: str, rve: Rve, output_file: str = "mesh_reqtri.mesh"
) -> None:
    """
    Prepares the mesh file for periodic remeshing

    :param mesh_file: mesh file to remesh
    :param rve: RVE for periodicity
    :param output_file: output file (must be .mesh)
    """
    gmsh.initialize()
    gmsh.open(mesh_file)

    gmsh.model.mesh.createTopology()  # Creates boundary entities

    # Encapsulate this block into a private function
    (
        _,
        surface_triangles_tags,
        surface_triangles_nodes_tags,
    ) = gmsh.model.mesh.getElements(_BOUNDARY_DIM)
    surface_triangles_nodes_tags: list[int] = list(surface_triangles_nodes_tags[0])
    surface_triangles_tags: list[int] = list(surface_triangles_tags[0])
    #####

    surface_triangles_nodes_coords = _get_surface_nodes_coords(
        surface_triangles_nodes_tags
    )

    surface_triangles = _build_surface_triangles(
        surface_triangles_tags, surface_triangles_nodes_coords
    )
    boundary_triangles_tags = _extract_boundary_triangles_tags(surface_triangles, rve)

    n_boundary_triangles = len(boundary_triangles_tags)

    # write files (there must be a more efficient way)

    gmsh.write(mesh_file[:-5] + "_triangles.mesh")
    gmsh.finalize()

    with open(mesh_file[:-5] + "_triangles.mesh", "r") as file:
        lines = file.readlines()

    with open(output_file, "w+") as outfile:
        outfile.writelines(lines[:-1])
        outfile.write("RequiredTriangles\n")
        outfile.write(str(n_boundary_triangles) + "\n")
        for tag in boundary_triangles_tags:
            outfile.write(str(tag) + "\n")
        outfile.write("End\n")


def _build_surface_triangles(
    surface_triangles_tags: list[int], surface_triangles_nodes_coords: np.ndarray
) -> list[Triangle]:
    n_triangles = len(surface_triangles_tags)
    surface_triangles = [
        Triangle(
            surface_triangles_nodes_coords[_MESH_DIM * i],
            surface_triangles_nodes_coords[_MESH_DIM * i + 1],
            surface_triangles_nodes_coords[_MESH_DIM * i + 2],
            surface_triangles_tags[i],
        )
        for i in range(n_triangles)
    ]
    return surface_triangles


def _compute_normal(triangle: Triangle) -> np.ndarray:
    """Computes normal vector of a triangle in mesh, defined by the coordinates of its three vertices"""
    u = triangle.node2 - triangle.node1
    v = triangle.node3 - triangle.node1
    n = np.cross(u, v)
    normal = n / np.linalg.norm(n)
    return normal


def _get_surface_nodes_coords(surface_nodes: list[int]) -> np.ndarray:
    n_nodes = len(surface_nodes)
    surface_nodes_coords = np.zeros((n_nodes, _MESH_DIM))
    for i in range(n_nodes):
        surface_nodes_coords[i] = gmsh.model.mesh.getNode(surface_nodes[i])[0]
    return surface_nodes_coords


def _is_triangle_on_boundary(triangle: Triangle, rve: Rve) -> bool:
    """Determines whether a triangle (defined by its 3 nodes) is on the boundary of a parallelepipedic rve"""

    xmin_normal = np.array([-1.0, 0.0, 0.0])
    ymin_normal = np.array([0.0, -1.0, 0.0])
    zmin_normal = np.array([0.0, 0.0, -1.0])
    xmax_normal = -xmin_normal
    ymax_normal = -ymin_normal
    zmax_normal = -zmin_normal

    # there must be a better way:

    bool_xmin = (
        m.isclose(triangle.node1[0], rve.x_min)
        and m.isclose(triangle.node2[0], rve.x_min)
        and m.isclose(triangle.node3[0], rve.x_min)
        and np.allclose(_compute_normal(triangle), xmin_normal)
    )

    bool_xmax = (
        m.isclose(triangle.node1[0], rve.x_max)
        and m.isclose(triangle.node2[0], rve.x_max)
        and m.isclose(triangle.node3[0], rve.x_max)
        and np.allclose(_compute_normal(triangle), xmax_normal)
    )

    bool_ymin = (
        m.isclose(triangle.node1[1], rve.y_min)
        and m.isclose(triangle.node2[1], rve.y_min)
        and m.isclose(triangle.node3[1], rve.y_min)
        and np.allclose(_compute_normal(triangle), ymin_normal)
    )

    bool_ymax = (
        m.isclose(triangle.node1[1], rve.y_max)
        and m.isclose(triangle.node2[1], rve.y_max)
        and m.isclose(triangle.node3[1], rve.y_max)
        and np.allclose(_compute_normal(triangle), ymax_normal)
    )

    bool_zmin = (
        m.isclose(triangle.node1[2], rve.z_min)
        and m.isclose(triangle.node2[2], rve.z_min)
        and m.isclose(triangle.node3[2], rve.z_min)
        and np.allclose(_compute_normal(triangle), zmin_normal)
    )

    bool_zmax = (
        m.isclose(triangle.node1[2], rve.z_max)
        and m.isclose(triangle.node2[2], rve.z_max)
        and m.isclose(triangle.node3[2], rve.z_max)
        and np.allclose(_compute_normal(triangle), zmax_normal)
    )

    return bool_xmin or bool_xmax or bool_ymin or bool_ymax or bool_zmin or bool_zmax


def _extract_boundary_triangles_tags(
    surface_triangles: list[Triangle], rve: Rve
) -> list[int]:
    n_tetrahedra = len(gmsh.model.mesh.getElements(_MESH_DIM)[1][0])
    boundary_triangles_tags = [
        int(triangle.tag - n_tetrahedra) # TODO add explicit comment
        for triangle in surface_triangles
        if _is_triangle_on_boundary(triangle, rve)
    ]

    return boundary_triangles_tags
