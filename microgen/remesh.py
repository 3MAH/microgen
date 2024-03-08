import os
import subprocess
from tempfile import NamedTemporaryFile
from typing import List, Optional, Union

import pyvista as pv

from microgen import BoxMesh, is_periodic


class InputBoxMeshNotPeriodicError(Exception):
    """Raised when the input BoxMesh is not periodic"""


class OutputMeshNotPeriodicError(Exception):
    """Raised when output mesh of remesh_keeping_periodicity_for_fem is not periodic"""


def remesh_keeping_periodicity_for_fem(
    input_mesh: BoxMesh,
    output_mesh_file: str,
    mesh_version: int = 2,
    dimension: int = 3,
    tol: float = 1e-8,
    hausd: Optional[float] = None,
    hgrad: Optional[float] = None,
    hmax: Optional[float] = None,
    hmin: Optional[float] = None,
    hsiz: Optional[float] = None,
) -> None:
    """
    Remeshes a mesh derived from a BoxMesh object using mmg while keeping periodicity

    :param input_mesh: BoxMesh to be remeshed
    :param output_mesh_file: output file (must be .mesh)
    :param mesh_version: mesh file version (default: 2)
    :param dimension: mesh dimension (default: 3)
    :param tol: tolerance for periodicity check

    The following parameters are used to control mmg remeshing, see here for more info : https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options
    :param hausd: Maximal Hausdorff distance for the boundaries approximation
    :param hgrad: Gradation value, ie ratio between lengths of adjacent mesh edges
    :param hmax: Maximal edge size
    :param hmin: Minimal edge size
    :param hsiz: Build a constant size map of size hsiz
    """
    with NamedTemporaryFile(
        suffix=".mesh", delete=False
    ) as boundary_triangles_file, NamedTemporaryFile(
        suffix=".mesh", delete=False
    ) as intermediate_mesh_file, NamedTemporaryFile(
        suffix=".mesh", delete=False
    ) as raw_output_mesh_file:
        _check_mesh_periodicity(input_mesh, tol, dimension)
        _generate_mesh_with_required_triangles(input_mesh, boundary_triangles_file.name)
        _preremesh_mmg(boundary_triangles_file.name, intermediate_mesh_file.name)
        _remesh_mmg(
            input_mesh_file=intermediate_mesh_file.name,
            output_mesh_file=raw_output_mesh_file.name,
            hausd=hausd,
            hgrad=hgrad,
            hmax=hmax,
            hmin=hmin,
            hsiz=hsiz,
        )
        os.remove(
            boundary_triangles_file.name
        )  # to solve compatibility issues of NamedTemporaryFiles with Windows
        os.remove(intermediate_mesh_file.name)
    os.remove(
        raw_output_mesh_file.name.replace(".mesh", ".sol")
    )  # Remove unused .sol file created by mmg
    _remove_unnecessary_fields_from_mesh_file(
        raw_output_mesh_file.name, output_mesh_file, mesh_version, dimension
    )
    _check_mesh_periodicity(output_mesh_file, tol, dimension)


def _check_mesh_periodicity(
    mesh: Union[str, BoxMesh], tol: float, dimension: int
) -> None:
    if isinstance(mesh, BoxMesh):
        nodes_coord = mesh.to_pyvista().points
    elif isinstance(mesh, str):
        pv_mesh = pv.read(mesh)
        nodes_coord = pv_mesh.points
    else:
        raise TypeError(
            "mesh parameter is neither a BoxMesh object or a mesh file string"
        )
    periodicity = is_periodic(nodes_coords=nodes_coord, tol=tol, dim=dimension)

    if not periodicity:
        raise InputBoxMeshNotPeriodicError("Input mesh is not periodic")


def _generate_mesh_with_required_triangles(
    input_mesh: BoxMesh, mesh_including_required_triangles: str = "merged_reqtri.mesh"
) -> None:
    with NamedTemporaryFile(suffix=".mesh", delete=True) as mesh_file:
        _generate_mesh_with_boundary_triangles(input_mesh, mesh_file.name)
        _add_required_triangles_to_mesh_file(
            input_mesh, mesh_file.name, mesh_including_required_triangles
        )


def _generate_mesh_with_boundary_triangles(
    input_mesh: BoxMesh, output_mesh: str = "merged.mesh"
) -> None:
    pyvista_mesh = input_mesh.to_pyvista()
    mesh_boundary, _ = input_mesh.boundary_elements(input_mesh.rve)
    merged_mesh = pyvista_mesh.merge(mesh_boundary)
    pv.save_meshio(output_mesh, merged_mesh)


def _get_number_of_boundary_triangles_from_boxmesh(input_mesh: BoxMesh) -> int:
    mesh_boundary, _ = input_mesh.boundary_elements(input_mesh.rve)
    n_boundary_triangles = mesh_boundary.n_cells

    return n_boundary_triangles


def _add_required_triangles_to_mesh_file(
    input_mesh: BoxMesh, input_mesh_file: str, output_mesh_file: str
) -> None:
    n_required_triangles = _get_number_of_boundary_triangles_from_boxmesh(input_mesh)
    with open(input_mesh_file) as input_file:
        lines = input_file.readlines()[:-1]  # remove last line End

    with open(output_mesh_file, "w+") as output_file:
        output_file.writelines(lines)
        output_file.write("RequiredTriangles\n")
        output_file.write(str(n_required_triangles) + "\n")
        for i in range(n_required_triangles):
            output_file.write(str(i + 1) + "\n")
        output_file.write("End\n")


def _remove_unnecessary_fields_from_mesh_file(
    input_mesh_file: str, output_mesh_file: str, mesh_version: int, dimension: int
) -> None:
    with open(input_mesh_file) as input_file:
        lines = input_file.readlines()

    write_bool = True
    with open(output_mesh_file, "w+") as output_file:
        output_file.write("MeshVersionFormatted " + str(mesh_version) + "\n\n")
        output_file.write("Dimension " + str(dimension) + "\n\n")
        for line in lines:
            if not _only_numbers_in_line(line.strip().split(" ")):
                if "Vertices" == line.strip() or "Tetrahedra" == line.strip():
                    write_bool = True
                else:
                    write_bool = False
            if write_bool:
                output_file.write(line)
        output_file.write("End\n")


def _only_numbers_in_line(str_list: List[str]) -> bool:
    return all(not flag.isalpha() for flag in str_list)


def _preremesh_mmg(
    input_mesh_file: str,
    output_mesh_file: str,
) -> None:
    mmg_system_call = [
        "mmg3d_O3",
        "-in",
        input_mesh_file,
        "-out",
        output_mesh_file,
        "-nofem",
    ]
    subprocess.call(mmg_system_call)


def _remesh_mmg(
    input_mesh_file: str,
    output_mesh_file: str,
    hausd: Optional[float] = None,
    hgrad: Optional[float] = None,
    hmax: Optional[float] = None,
    hmin: Optional[float] = None,
    hsiz: Optional[float] = None,
) -> None:
    mmg_system_call = [
        "mmg3d_O3",
        "-in",
        input_mesh_file,
        "-out",
        output_mesh_file,
        "-ls",
        "-nr",
    ]
    if hausd:
        mmg_system_call.extend(["-hausd", str(hausd)])
    if hgrad:
        mmg_system_call.extend(["-hgrad", str(hgrad)])
    if hmax:
        mmg_system_call.extend(["-hmax", str(hmax)])
    if hmin:
        mmg_system_call.extend(["-hmin", str(hmin)])
    if hsiz:
        mmg_system_call.extend(["-hsiz", str(hsiz)])
    subprocess.call(mmg_system_call)
