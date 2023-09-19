import os
from tempfile import NamedTemporaryFile

from microgen import BoxMesh, Mmg
import meshio

def remesh_keeping_periodicity_for_fem(
    input_mesh: BoxMesh,
    output_mesh_file: str,
    hausd: float = None,
    hgrad: float = None,
    hmax: float = None,
    hmin: float = None,
    hsiz: float = None,
) -> None:
    """
    Remeshes a mesh (.mesh file format) derived from a BoxMesh using mmg while keeping periodicity

    :param input_mesh: BoxMesh to be remeshed
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
        _generate_mesh_with_required_triangles(input_mesh, boundary_triangles_file.name)
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
    os.remove(output_mesh_file.replace(".mesh",".sol")) # Remove unused .sol file created by mmg


def _generate_mesh_with_required_triangles(input_mesh : BoxMesh, mesh_including_required_triangles : str = 'merged_reqtri.mesh') -> None:
    with NamedTemporaryFile(suffix='.vtk', delete=False) as vtk_file, NamedTemporaryFile(suffix='.mesh', delete=False) as mesh_file:
        _generate_vtk_with_boundary_triangles(input_mesh, vtk_file.name)
        _convert_vtk_to_mesh(vtk_file.name, mesh_file.name)
        _add_required_triangles_to_mesh_file(input_mesh, mesh_file.name, mesh_including_required_triangles)


def _generate_vtk_with_boundary_triangles(input_mesh : BoxMesh, output_mesh : str = 'merged.vtk') -> None:
    pyvista_mesh = input_mesh.to_pyvista()
    mesh_boundary, _ = input_mesh.boundary_elements(input_mesh.rve)
    merged_mesh = pyvista_mesh.merge(mesh_boundary)
    merged_mesh.save(output_mesh, binary=False)

def _convert_vtk_to_mesh(input_vtk_file : str, output_mesh_file : str) -> None:
    meshio_mesh = meshio.read(input_vtk_file)
    meshio_mesh.write(output_mesh_file)


def _get_number_of_boundary_triangles_from_boxmesh(input_mesh : BoxMesh) -> int:
    mesh_boundary, _ = input_mesh.boundary_elements(input_mesh.rve)
    n_boundary_triangles = mesh_boundary.n_faces

    return n_boundary_triangles

def _add_required_triangles_to_mesh_file(input_mesh : BoxMesh, input_mesh_file : str, output_mesh_file : str) -> None:
    n_required_triangles = _get_number_of_boundary_triangles_from_boxmesh(input_mesh)
    with open(input_mesh_file, 'r') as input_file:
        lines = input_file.readlines()[:-1] #remove last line End

    with open(output_mesh_file, 'w+') as output_file:
        output_file.writelines(lines)
        output_file.write("RequiredTriangles\n")
        output_file.write(str(n_required_triangles) + "\n")
        for i in range(n_required_triangles):
            output_file.write(str(i+1) + "\n")
        output_file.write("End\n")
