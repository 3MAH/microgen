import math as m
import os
from tempfile import NamedTemporaryFile
from typing import NamedTuple

import gmsh
import numpy as np
import numpy.typing as npt

from microgen import Mmg, Rve


class Triangle(NamedTuple):
    node1: npt.NDArray[np.float_]
    node2: npt.NDArray[np.float_]
    node3: npt.NDArray[np.float_]
    tag: int


_MESH_DIM = 3
_BOUNDARY_DIM = 2


def remesh_keeping_periodicity_for_fem(
    input_mesh_file: str,
    rve: Rve,
    output_mesh_file: str,
    hausd: float = None,
    hgrad: float = None,
    hmax: float = None,
    hmin: float = None,
    hsiz: float = None,
) -> None:
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
    with NamedTemporaryFile(
        suffix=".mesh", delete=False
    ) as boundary_triangles_file:
        _identify_boundary_triangles_from_mesh_file(
            input_mesh_file, rve, boundary_triangles_file.name
        )
        Mmg.mmg3d(
            input=boundary_triangles_file.name,
            output=output_mesh_file,
            hausd=hausd,
            hgrad=hgrad,
            hmax=hmax,
            hmin=hmin,
            hsiz=hsiz,
        )
        os.remove(
            boundary_triangles_file.name
        )  # to solve compatibility issues of NamedTemporaryFiles with Windows


def _identify_boundary_triangles_from_mesh_file(
    input_mesh_file: str, rve: Rve, output_mesh_file: str
) -> None:
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
    boundary_triangles_tags: list[int] = _extract_boundary_triangles_tags(
        surface_triangles, rve
    )

    _write_output(boundary_triangles_tags, output_mesh_file)


def _get_surface_triangles() -> tuple[list[int], list[int]]:
    gmsh.model.mesh.createTopology()  # Creates boundary entities
    (
        _,
        surface_triangles_tags,
        surface_triangles_nodes_tags,
    ) = gmsh.model.mesh.getElements(_BOUNDARY_DIM)
    surface_triangles_nodes_tags: list[int] = list(
        surface_triangles_nodes_tags[0]
    )
    surface_triangles_tags: list[int] = list(surface_triangles_tags[0])

    return surface_triangles_tags, surface_triangles_nodes_tags


def _build_surface_triangles() -> list[Triangle]:
    (
        surface_triangles_tags,
        surface_triangles_nodes_tags,
    ) = _get_surface_triangles()
    surface_triangles_nodes_coords = _get_surface_nodes_coords(
        surface_triangles_nodes_tags
    )

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


def _get_surface_nodes_coords(
    surface_nodes: list[int],
) -> npt.NDArray[np.float_]:
    n_nodes = len(surface_nodes)
    surface_nodes_coords = np.zeros((n_nodes, _MESH_DIM))
    for i in range(n_nodes):
        surface_nodes_coords[i] = gmsh.model.mesh.getNode(surface_nodes[i])[0]
    return surface_nodes_coords


def _is_triangle_on_boundary(triangle: Triangle, rve: Rve) -> bool:
    """Determines whether a triangle (defined by its 3 nodes) is on the boundary of a parallelepipedic rve"""

    rve_boundaries = [
        (rve.x_min, rve.x_max),
        (rve.y_min, rve.y_max),
        (rve.z_min, rve.z_max),
    ]
    triangle_nodes = triangle.node1, triangle.node2, triangle.node3
    for i, rve_axis_min_max in enumerate(rve_boundaries):
        for rve_axis_boundary in rve_axis_min_max:
            is_triangle_on_boundary = all(
                m.isclose(node[i], rve_axis_boundary) for node in triangle_nodes
            )
            if is_triangle_on_boundary:
                return True
    return False


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


def _write_output(
    boundary_triangles_tags: list[int], output_file: str = "mesh_reqtri.mesh"
) -> None:
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
