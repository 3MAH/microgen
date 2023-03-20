import gmsh
from microgen import Rve, Mmg
import numpy as np
import math as m
from tempfile import NamedTemporaryFile

from typing import NamedTuple
import numpy.typing as npt


class Triangle(NamedTuple):
    node1: npt.NDArray[np.float_]
    node2: npt.NDArray[np.float_]
    node3: npt.NDArray[np.float_]
    tag: int

_MESH_DIM = 3
_BOUNDARY_DIM = 2


def remesh_keeping_periodicity(input_mesh_file: str, rve: Rve,
                               output_mesh_file: str, hausd: float = None,
                               hgrad: float = None,
                               hmax: float = None,
                               hmin: float = None,
                               hsiz: float = None) -> None:
    """
    Remeshes a mesh (.mesh file format) using mmg while keeping periodicity

    :param input_mesh_file: mesh file to remesh (must be .mesh)
    :param rve: Representative Volume Element for periodicity
    :param output_mesh_file: output file (must be .mesh)

    The following parameters are used to control mmg remeshing, see here for more info : https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options
    :param hausd: Maximal Hausdorff distance for the boundaries approximation
    :param hgrad: Gradation value, ie ratio between lengths of adjacent mesh edges
    :param hmax: Maximal edge size
    :param hmin: Minimal edge size
    :param hsiz: Build a constant size map of size hsiz
    """
    with NamedTemporaryFile(suffix='.mesh') as boundary_triangles_file:
        _identify_boundary_triangles_from_mesh_file(input_mesh_file, rve, boundary_triangles_file.name)
        Mmg.mmg3d(input=boundary_triangles_file.name, output=output_mesh_file, hausd=hausd, hgrad=hgrad, hmax=hmax,
                  hmin=hmin, hsiz=hsiz)


def _identify_boundary_triangles_from_mesh_file(input_mesh_file: str, rve: Rve,
                                                output_mesh_file: str) -> None:
    """
    Prepares the mesh file for periodic remeshing by tagging boundary triangles

    This function adds a "RequiredTriangles" field
    in the ".mesh" file that stores the boundary triangles tags.
    mmg recognizes this field as triangles it must not remesh

    :param input_mesh_file: mesh file to remesh (must be .mesh)
    :param rve: Representative Volume Element for periodicity
    :param output_mesh_file: output file (must be .mesh)
    """
    gmsh.initialize()
    gmsh.open(input_mesh_file)

    surface_triangles: list[Triangle] = _build_surface_triangles()
    boundary_triangles_tags: list[int] = _extract_boundary_triangles_tags(surface_triangles, rve)

    _write_output(boundary_triangles_tags, output_mesh_file)


def _get_surface_triangles() -> tuple[list[int], list[int]]:
    gmsh.model.mesh.createTopology()  # Creates boundary entities
    (
        _,
        surface_triangles_tags,
        surface_triangles_nodes_tags,
    ) = gmsh.model.mesh.getElements(_BOUNDARY_DIM)
    surface_triangles_nodes_tags: list[int] = list(surface_triangles_nodes_tags[0])
    surface_triangles_tags: list[int] = list(surface_triangles_tags[0])

    return surface_triangles_tags, surface_triangles_nodes_tags


def _build_surface_triangles() -> list[Triangle]:
    surface_triangles_tags, surface_triangles_nodes_tags = _get_surface_triangles()
    surface_triangles_nodes_coords = _get_surface_nodes_coords(
        surface_triangles_nodes_tags)

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


def _get_surface_nodes_coords(surface_nodes: list[int]) -> np.ndarray:
    n_nodes = len(surface_nodes)
    surface_nodes_coords = np.zeros((n_nodes, _MESH_DIM))
    for i in range(n_nodes):
        surface_nodes_coords[i] = gmsh.model.mesh.getNode(surface_nodes[i])[0]
    return surface_nodes_coords


def _is_triangle_on_boundary(triangle: Triangle, rve: Rve) -> bool:
    """Determines whether a triangle (defined by its 3 nodes) is on the boundary of a parallelepipedic rve"""

    # there must be a better way:

    bool_xmin = (
            m.isclose(triangle.node1[0], rve.x_min)
            and m.isclose(triangle.node2[0], rve.x_min)
            and m.isclose(triangle.node3[0], rve.x_min)
    )

    bool_xmax = (
            m.isclose(triangle.node1[0], rve.x_max)
            and m.isclose(triangle.node2[0], rve.x_max)
            and m.isclose(triangle.node3[0], rve.x_max)
    )

    bool_ymin = (
            m.isclose(triangle.node1[1], rve.y_min)
            and m.isclose(triangle.node2[1], rve.y_min)
            and m.isclose(triangle.node3[1], rve.y_min)
    )

    bool_ymax = (
            m.isclose(triangle.node1[1], rve.y_max)
            and m.isclose(triangle.node2[1], rve.y_max)
            and m.isclose(triangle.node3[1], rve.y_max)
    )

    bool_zmin = (
            m.isclose(triangle.node1[2], rve.z_min)
            and m.isclose(triangle.node2[2], rve.z_min)
            and m.isclose(triangle.node3[2], rve.z_min)
    )

    bool_zmax = (
            m.isclose(triangle.node1[2], rve.z_max)
            and m.isclose(triangle.node2[2], rve.z_max)
            and m.isclose(triangle.node3[2], rve.z_max)
    )

    return bool_xmin or bool_xmax or bool_ymin or bool_ymax or bool_zmin or bool_zmax


def _extract_boundary_triangles_tags(
        surface_triangles: list[Triangle], rve: Rve
) -> list[int]:
    n_tetrahedra = len(gmsh.model.mesh.getElements(_MESH_DIM)[1][0])
    boundary_triangles_tags = [
        int(triangle.tag - n_tetrahedra)
        # to make sure gmsh and mmg count elements the same way (ie: do not confuse tetrahedra and triangles)
        for triangle in surface_triangles
        if _is_triangle_on_boundary(triangle, rve)
    ]

    return boundary_triangles_tags


def _write_output(boundary_triangles_tags: list[int], output_file: str = "mesh_reqtri.mesh") -> None:
    n_boundary_triangles = len(boundary_triangles_tags)
    gmsh.write(output_file)
    gmsh.finalize()

    with open(output_file, "r") as file:
        lines = file.readlines()[:-1]  # Delete last line End

    with open(output_file, "w") as outfile:
        outfile.writelines(lines)
        outfile.write("RequiredTriangles\n")
        outfile.write(str(n_boundary_triangles) + "\n")
        for tag in boundary_triangles_tags:
            outfile.write(str(tag) + "\n")
        outfile.write("End\n")
