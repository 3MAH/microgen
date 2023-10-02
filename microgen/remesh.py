import subprocess

from microgen import BoxMesh, Rve
from typing import Optional
import meshio

import pyvista as pv


def remesh_keeping_periodicity_for_fem(
        input_mesh: pv.UnstructuredGrid,
        rve: Optional[Rve],
        output_mesh_file: str,
        boundary_triangles_file = 'merged_reqtri.mesh',
        raw_output_mesh_file = 'raw_output_mesh.mesh',
        mesh_version : int = 2,
        dimension : int = 3,
        hausd: float = None,
        hgrad: float = None,
        hmax: float = None,
        hmin: float = None,
        hsiz: float = None,
) -> None:
    """
    Remeshes a mesh (.mesh file format) derived from a BoxMesh using mmg while keeping periodicity

    :param input_mesh: BoxMesh to be remeshed
    :param rve: associated Representative Volume Element (bounding box of mesh)
    :param output_mesh_file: output file (must be .mesh)

    :param boundary_triangles_file: .mesh file containing the required triangles field
    :param raw_output_mesh_file: optimized .mesh file with unused fields blocking meshio conversion
    :param mesh_version: .mesh file version (only version 2 supported for now)
    :param dimension: mesh dimension (only 3D tet meshes supported for now)

    The following parameters are used to control mmg remeshing, see here for more info : https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options
    :param hausd: Maximal Hausdorff distance for the boundaries approximation
    :param hgrad: Gradation value, ie ratio between lengths of adjacent mesh edges
    :param hmax: Maximal edge size
    :param hmin: Minimal edge size
    :param hsiz: Build a constant size map of size hsiz
    """
    _generate_mesh_with_required_triangles(input_mesh, rve, boundary_triangles_file)
    _remesh_mmg(
        input_mesh_file=boundary_triangles_file,
        output_mesh_file=raw_output_mesh_file,
        hausd=hausd,
        hgrad=hgrad,
        hmax=hmax,
        hmin=hmin,
        hsiz=hsiz,
    )
    _remove_unnecessary_fields_from_mesh_file(raw_output_mesh_file, output_mesh_file, mesh_version, dimension)


def _generate_mesh_with_required_triangles(input_mesh: pv.UnstructuredGrid, rve: Optional[Rve] = None, mesh_including_required_triangles : str = 'mesh_reqtri.mesh') -> None:
    input_boxmesh = BoxMesh.from_pyvista(input_mesh)
    if rve is not None:
        input_boxmesh.rve = rve
    mesh_boundary, _ = input_boxmesh.boundary_elements(input_boxmesh.rve)
    mesh_boundary.save('boundary.vtk')
    pv.save_meshio('boundary.mesh', mesh_boundary)
    merged_mesh = input_mesh.merge(mesh_boundary)
    merged_mesh.save('input_with_triangles.vtk')
    pv.save_meshio('input_with_triangles.mesh', merged_mesh)

    n_boundary_triangles = mesh_boundary.n_faces

    with open('input_with_triangles.mesh', 'r') as merged_mesh_with_boundary_file:
        lines = merged_mesh_with_boundary_file.readlines()[:-1] # remove last line End

    with open(mesh_including_required_triangles, 'w') as mesh_with_required_triangles_file:
        mesh_with_required_triangles_file.writelines(lines)
        mesh_with_required_triangles_file.write("RequiredTriangles\n")
        mesh_with_required_triangles_file.write(str(n_boundary_triangles) + "\n")
        for i in range(n_boundary_triangles):
            mesh_with_required_triangles_file.write(str(i + 1) + "\n")
        mesh_with_required_triangles_file.write("End\n")



def _remove_unnecessary_fields_from_mesh_file(input_mesh_file: str, output_mesh_file: str, mesh_version: int, dimension: int) -> None:
    with open(input_mesh_file, 'r') as input_file:
        lines = input_file.readlines()

    write_bool = True
    with open(output_mesh_file, 'w+') as output_file:
        output_file.write('MeshVersionFormatted ' + str(mesh_version) + '\n\n')
        output_file.write('Dimension ' + str(dimension) + '\n\n')
        for line in lines:
            if not _only_numbers_in_line(line.strip().split(' ')):
                if 'Vertices' == line.strip() or 'Tetrahedra' == line.strip():
                    write_bool = True
                else:
                    write_bool = False
            if write_bool:
                output_file.write(line)
        output_file.write('End\n')


def _only_numbers_in_line(str_list: list[str]) -> bool:
    return all(not flag.isalpha() for flag in str_list)

def _remesh_mmg(
    input_mesh_file: str,
    output_mesh_file: str,
    hausd : float = None,
    hgrad : float = None,
    hmax : float = None,
    hmin : float = None,
    hsiz : float = None,
) -> None:
    mmg_system_call = [
        "mmg3d_O3", "-in",
        input_mesh_file,
        "-out",
        output_mesh_file,
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