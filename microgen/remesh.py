"""Remesh a mesh using mmg while keeping boundaries untouched."""

from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import overload

import numpy as np
import numpy.typing as npt
import pyvista as pv

from microgen import BoxMesh, Mmg, is_periodic
from microgen.external import MmgError

# mmg3d's iso-zero adaptation (`-ls 0`) is non-deterministic on the metric .sol
# left by the first pass and can produce a non-periodic mesh, or crash, on
# dense / high-frequency input. Empirically the next attempt almost always
# converges, so retry the whole pipeline up to this many times before giving
# up.
_MMG_MAX_ATTEMPTS = 5


class InputMeshNotPeriodicError(Exception):
    """Raised when the input mesh is not periodic.

    Triggered when ``periodic=True`` is passed to
    :func:`remesh_keeping_boundaries_for_fem` but the supplied input mesh
    fails the periodicity check.
    """


class OutputMeshNotPeriodicError(Exception):
    """Raised when the remeshed output is not periodic.

    Triggered when ``periodic=True`` is passed to
    :func:`remesh_keeping_boundaries_for_fem` but mmg's output fails the
    periodicity check.
    """


@overload
def remesh_keeping_boundaries_for_fem(
    input_mesh: BoxMesh,
    *,
    periodic: bool = True,
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
def remesh_keeping_boundaries_for_fem(
    input_mesh: pv.UnstructuredGrid,
    *,
    periodic: bool = True,
    mesh_version: int = 2,
    dimension: int = 3,
    tol: float = 1e-8,
    hausd: float | None = None,
    hgrad: float | None = None,
    hmax: float | None = None,
    hmin: float | None = None,
    hsiz: float | None = None,
) -> pv.UnstructuredGrid: ...


def remesh_keeping_boundaries_for_fem(  # noqa: PLR0913
    input_mesh: BoxMesh | pv.UnstructuredGrid,
    *,
    periodic: bool = True,
    mesh_version: int = 2,
    dimension: int = 3,
    tol: float = 1e-8,
    hausd: float | None = None,
    hgrad: float | None = None,
    hmax: float | None = None,
    hmin: float | None = None,
    hsiz: float | None = None,
) -> BoxMesh | pv.UnstructuredGrid:
    """Remesh a mesh using mmg while keeping boundary elements untouched.

    See https://www.mmgtools.org/mmg-remesher-try-mmg/mmg-remesher-options
    for the full reference on each ``h*`` parameter.

    :param input_mesh: BoxMesh or pv.UnstructuredGrid mesh to be remeshed
    :param periodic: whether the mesh is periodic and must stay periodic
    :param mesh_version: mesh file version
    :param dimension: mesh dimension
    :param tol: tolerance for the periodicity check
    :param hausd: maximal Hausdorff distance for the boundaries approximation
    :param hgrad: gradation value, i.e. ratio between lengths of adjacent
        mesh edges
    :param hmax: maximal edge size
    :param hmin: minimal edge size
    :param hsiz: build a constant size map of size ``hsiz``
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

    if periodic and not is_periodic(nodes_coords, tol):
        err_msg = "Input mesh is not periodic"
        raise InputMeshNotPeriodicError(err_msg)

    with NamedTemporaryFile(suffix=".mesh", delete=False) as boundary_triangles_file:
        pass

    _generate_mesh_with_required_triangles(
        input_box_mesh,
        boundary_triangles_file.name,
    )

    last_error: Exception | None = None
    output_mesh: pv.UnstructuredGrid | None = None
    for _ in range(_MMG_MAX_ATTEMPTS):
        try:
            output_mesh = _run_mmg_pipeline(
                boundary_triangles_file.name,
                mesh_version=mesh_version,
                dimension=dimension,
                hausd=hausd,
                hgrad=hgrad,
                hmax=hmax,
                hmin=hmin,
                hsiz=hsiz,
            )
        except MmgError as err:
            last_error = err
            continue

        if periodic:
            _snap_to_expected_bounds(output_mesh, nodes_coords, tol)
            if not is_periodic(output_mesh.points, tol):
                last_error = OutputMeshNotPeriodicError(
                    "Something went wrong: output mesh is not periodic",
                )
                output_mesh = None
                continue
        break

    Path(boundary_triangles_file.name).unlink(missing_ok=True)

    if output_mesh is None:
        assert last_error is not None
        raise last_error

    if isinstance(input_mesh, BoxMesh):
        return BoxMesh.from_pyvista(output_mesh)
    return output_mesh


def _run_mmg_pipeline(  # noqa: PLR0913
    boundary_triangles_file: str,
    *,
    mesh_version: int,
    dimension: int,
    hausd: float | None,
    hgrad: float | None,
    hmax: float | None,
    hmin: float | None,
    hsiz: float | None,
) -> pv.UnstructuredGrid:
    """Run the two mmg3d passes and return the parsed output mesh.

    Owns its own temp files so a caller retrying after a failure starts from
    a clean state.
    """
    with (
        NamedTemporaryFile(suffix=".mesh", delete=False) as premeshed_mesh_file,
        NamedTemporaryFile(suffix=".mesh", delete=False) as raw_output_mesh_file,
        NamedTemporaryFile(suffix=".mesh", delete=False) as output_mesh_file,
    ):
        pass
    try:
        Mmg.mmg3d(
            input=boundary_triangles_file,
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
        return pv.UnstructuredGrid(output_mesh_file.name)
    finally:
        for path in (
            premeshed_mesh_file.name,
            premeshed_mesh_file.name.replace(".mesh", ".sol"),
            raw_output_mesh_file.name,
            raw_output_mesh_file.name.replace(".mesh", ".sol"),
            output_mesh_file.name,
        ):
            Path(path).unlink(missing_ok=True)


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
        output_file.write(f"{n_required_triangles}\n")
        output_file.writelines(f"{i + 1}\n" for i in range(n_required_triangles))
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
        output_file.write(f"MeshVersionFormatted {mesh_version}\n\n")
        output_file.write(f"Dimension {dimension}\n\n")
        for line in lines:
            if not _only_numbers_in_line(line.strip().split(" ")):
                write_bool = line.strip() in ("Vertices", "Tetrahedra")
            if write_bool:
                output_file.write(line)
        output_file.write("End\n")


def _only_numbers_in_line(line: list[str]) -> bool:
    return all(not flag.isalpha() for flag in line)


def _snap_to_expected_bounds(
    mesh: pv.UnstructuredGrid,
    reference_coords: npt.NDArray[np.float64],
    tol: float,
) -> None:
    """Reconcile the output bounding box with the input's.

    mmg sometimes writes coordinates that drift past the input bounds:

    * by floating-point noise (e.g. y = -1e-22 instead of 0). The bulk of the
      boundary is intact; snap coordinates within ``tol`` of an expected face
      back exactly onto it.
    * by a few microns (e.g. y = 1.00006 instead of 1.0). The boundary face
      itself still has the correct node count, so this stray vertex does not
      have a periodic twin. Pull it just inside the cell -- it then contributes
      no boundary node and doesn't anchor ``_get_bounding_box`` past the cell.

    Both forms confuse ``is_periodic`` even when the periodic boundary is
    otherwise intact.
    """
    expected_min = reference_coords.min(axis=0)
    expected_max = reference_coords.max(axis=0)
    interior_offset = tol * 10.0
    points = np.asarray(mesh.points, dtype=np.float64).copy()
    for axis in range(points.shape[1]):
        col = points[:, axis]
        snap_min = np.abs(col - expected_min[axis]) <= tol
        snap_max = np.abs(col - expected_max[axis]) <= tol
        col[snap_min] = expected_min[axis]
        col[snap_max] = expected_max[axis]
        max_drift = tol * 100.0
        below = (col < expected_min[axis]) & (col >= expected_min[axis] - max_drift)
        above = (col > expected_max[axis]) & (col <= expected_max[axis] + max_drift)
        col[below] = expected_min[axis] + interior_offset
        col[above] = expected_max[axis] - interior_offset
        points[:, axis] = col
    mesh.points = points
