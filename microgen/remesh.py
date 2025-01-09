"""Remesh a mesh using mmg while keeping periodicity."""

from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import overload

import pyvista as pv

from microgen import BoxMesh, Mmg, is_periodic


class InputMeshNotPeriodicError(Exception):
    """Raised when input mesh of remesh_keeping_periodicity_for_fem is not periodic."""


class OutputMeshNotPeriodicError(Exception):
    """Raised when output mesh of remesh_keeping_periodicity_for_fem is not periodic."""


@overload
def remesh_keeping_periodicity_for_fem(
    input_mesh: BoxMesh,
    mesh_version: int = 2,
    dimension: int = 3,
    tol: float = 1e-8,
    hausd: float | None = None,
    hgrad: float | None = None,
    hmax: float | None = None,
    hmin: float | None = None,
    hsiz: float | None = None,
) -> BoxMesh: ...


@overload
def remesh_keeping_periodicity_for_fem(
    input_mesh: pv.UnstructuredGrid,
    mesh_version: int = 2,
    dimension: int = 3,
    tol: float = 1e-8,
    hausd: float | None = None,
    hgrad: float | None = None,
    hmax: float | None = None,
    hmin: float | None = None,
    hsiz: float | None = None,
) -> pv.UnstructuredGrid: ...


def remesh_keeping_periodicity_for_fem(
    input_mesh: BoxMesh | pv.UnstructuredGrid,
    mesh_version: int = 2,
    dimension: int = 3,
    tol: float = 1e-8,
    hausd: float | None = None,
    hgrad: float | None = None,
    hmax: float | None = None,
    hmin: float | None = None,
    hsiz: float | None = None,
) -> BoxMesh | pv.UnstructuredGrid:
    """Remesh a mesh using mmg while keeping periodicity.

    :param input_mesh: BoxMesh or pv.UnstructuredGrid mesh to be remeshed
    :param mesh_version: mesh file version (default: 2)
    :param dimension: mesh dimension (default: 3)
    :param tol: tolerance for periodicity check

    The following parameters are used to control mmg remeshing, see here for more info :
    https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options

    :param hausd: Maximal Hausdorff distance for the boundaries approximation
    :param hgrad: Gradation value, ie ratio between lengths of adjacent mesh edges
    :param hmax: Maximal edge size
    :param hmin: Minimal edge size
    :param hsiz: Build a constant size map of size hsiz
    """
    if isinstance(input_mesh, pv.UnstructuredGrid):
        is_only_tetra = (
            len(input_mesh.cells_dict) == 1
            and pv.CellType.TETRA in input_mesh.cells_dict
        )
        if not is_only_tetra:
            input_mesh = input_mesh.triangulate()
        nodes_coords = input_mesh.points
        input_box_mesh = BoxMesh.from_pyvista(input_mesh)
    elif isinstance(input_mesh, BoxMesh):
        nodes_coords = input_mesh.to_pyvista().points
        input_box_mesh = input_mesh
    else:
        err_msg = "Input mesh must be either a BoxMesh or a pv.UnstructuredGrid"
        raise TypeError(err_msg)

    if not is_periodic(nodes_coords, tol, dimension):
        err_msg = "Input mesh is not periodic"
        raise InputMeshNotPeriodicError(err_msg)

    with (
        NamedTemporaryFile(suffix=".mesh", delete=False) as boundary_triangles_file,
        NamedTemporaryFile(suffix=".mesh", delete=False) as premeshed_mesh_file,
        NamedTemporaryFile(suffix=".mesh", delete=False) as raw_output_mesh_file,
        NamedTemporaryFile(suffix=".mesh", delete=False) as output_mesh_file,
    ):
        _generate_mesh_with_required_triangles(
            input_box_mesh,
            boundary_triangles_file.name,
        )
        Mmg.mmg3d(
            input=boundary_triangles_file.name,
            output=premeshed_mesh_file.name,
            nofem=True,
        )
        Mmg.mmg3d(
            input=premeshed_mesh_file.name,
            output=raw_output_mesh_file.name,
            hausd=hausd,
            hgrad=hgrad,
            hmax=hmax,
            hmin=hmin,
            hsiz=hsiz,
            ls=True,
            nr=True,
        )

    _remove_unnecessary_fields_from_mesh_file(
        raw_output_mesh_file.name,
        output_mesh_file.name,
        mesh_version,
        dimension,
    )

    output_mesh = pv.UnstructuredGrid(output_mesh_file.name)

    if not is_periodic(output_mesh.points, tol, dimension):
        err_msg = "Something went wrong: output mesh is not periodic"
        raise OutputMeshNotPeriodicError(err_msg)

    # Remove unused .sol files created by mmg
    # Solve compatibility issues of NamedTemporaryFiles with Windows
    trash_files_list = [
        boundary_triangles_file.name,
        premeshed_mesh_file.name,
        premeshed_mesh_file.name.replace(".mesh", ".sol"),
        raw_output_mesh_file.name,
        raw_output_mesh_file.name.replace(".mesh", ".sol"),
        output_mesh_file.name,
    ]
    for file in trash_files_list:
        Path(file).unlink()

    if isinstance(input_mesh, BoxMesh):
        return BoxMesh.from_pyvista(output_mesh)
    return output_mesh


def _generate_mesh_with_required_triangles(
    input_mesh: BoxMesh,
    mesh_including_required_triangles: str = "merged_reqtri.mesh",
) -> None:
    with NamedTemporaryFile(suffix=".mesh", delete=True) as mesh_file:
        _generate_mesh_with_boundary_triangles(input_mesh, mesh_file.name)
        _add_required_triangles_to_mesh_file(
            input_mesh,
            mesh_file.name,
            mesh_including_required_triangles,
        )


def _generate_mesh_with_boundary_triangles(
    input_mesh: BoxMesh,
    output_mesh: str = "merged.mesh",
) -> None:
    pyvista_mesh = input_mesh.to_pyvista()
    mesh_boundary, _ = input_mesh.boundary_elements(input_mesh.rve)
    merged_mesh = pyvista_mesh.merge(mesh_boundary)
    pv.save_meshio(output_mesh, merged_mesh)


def _get_number_of_boundary_triangles_from_boxmesh(input_mesh: BoxMesh) -> int:
    mesh_boundary, _ = input_mesh.boundary_elements(input_mesh.rve)
    return mesh_boundary.n_cells


def _add_required_triangles_to_mesh_file(
    input_mesh: BoxMesh,
    input_mesh_file: str,
    output_mesh_file: str,
) -> None:
    n_required_triangles = _get_number_of_boundary_triangles_from_boxmesh(input_mesh)
    with Path(input_mesh_file).open() as input_file:
        lines = input_file.readlines()[:-1]  # remove last line End

    with Path(output_mesh_file).open(mode="w+") as output_file:
        output_file.writelines(lines)
        output_file.write("RequiredTriangles\n")
        output_file.write(str(n_required_triangles) + "\n")
        for i in range(n_required_triangles):
            output_file.write(str(i + 1) + "\n")
        output_file.write("End\n")


def _remove_unnecessary_fields_from_mesh_file(
    input_mesh_file: str,
    output_mesh_file: str,
    mesh_version: int,
    dimension: int,
) -> None:
    with Path(input_mesh_file).open() as input_file:
        lines = input_file.readlines()

    write_bool = True
    with Path(output_mesh_file).open(mode="w+") as output_file:
        output_file.write(f"MeshVersionFormatted {mesh_version}" + "\n\n")
        output_file.write(f"Dimension {dimension}" + "\n\n")
        for line in lines:
            if not _only_numbers_in_line(line.strip().split(" ")):
                write_bool = line.strip() in ("Vertices", "Tetrahedra")
            if write_bool:
                output_file.write(line)
        output_file.write("End\n")


def _only_numbers_in_line(line: list[str]) -> bool:
    return all(not flag.isalpha() for flag in line)
