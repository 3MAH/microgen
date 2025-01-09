"""Tests for the meshPeriodic function."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import cadquery as cq
import numpy as np
import pytest
import pyvista as pv

from microgen import (
    Box,
    Cylinder,
    Phase,
    Rve,
    cutPhases,
    fuseShapes,
    is_periodic,
    meshPeriodic,
    periodic,
)

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/
# ruff: noqa: E501 line-too-long https://docs.astral.sh/ruff/rules/line-too-long/


@pytest.fixture()
def rve_unit() -> Rve:
    """Return a unit Rve."""
    return Rve(dim=1, center=(0.5, 0.5, 0.5))


@pytest.fixture()
def rve_double() -> Rve:
    """Return a double sized Rve."""
    return Rve(dim=2, center=(1.0, 1.0, 1.0))


@pytest.fixture()
def rve_double_centered() -> Rve:
    """Return a double sized Rve centered at the origin."""
    return Rve(dim=2, center=(0.0, 0.0, 0.0))


@pytest.fixture()
def tmp_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Return a temporary directory path."""
    tmp_dir_name = "test_tmp_dir"
    return tmp_path_factory.mktemp(tmp_dir_name)


@pytest.fixture()
def tmp_output_compound_filename(tmp_dir: Path) -> str:
    """Return a temporary compound filename."""
    return (tmp_dir / "compound.step").as_posix()


@pytest.fixture()
def tmp_output_vtk_filename(tmp_dir: Path) -> str:
    """Return a temporary VTK filename."""
    return (tmp_dir / "octettruss.vtk").as_posix()


def _generate_cqcompound_octettruss(rve: Rve) -> list[Phase]:
    """Generate a list of periodic octet truss cylinders."""
    # data to generate octet truss cylinders
    xc = np.array([0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0])
    yc = np.array([0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5])
    zc = np.array([0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5])
    psi = np.array(
        [45.0, -45.0, 0.0, 0.0, 90.0, 90.0, 0.0, 0.0, 90.0, 90.0, 45.0, -45.0],
    )
    theta = np.array(
        [0.0, 0.0, -90.0, -90.0, -90.0, -90.0, 90.0, 90.0, 90.0, 90.0, 0.0, 0.0],
    )
    phi = np.array(
        [0.0, 0.0, 45.0, -45.0, 45.0, -45.0, -45.0, 45.0, -45.0, 45.0, 0.0, 0.0],
    )
    height = np.full_like(xc, 2.0)
    radius = np.full_like(xc, 0.05)

    phases: list[Phase] = []

    for i in range(len(xc)):
        elem = Cylinder(
            center=(xc[i], yc[i], zc[i]),
            orientation=(psi[i], theta[i], phi[i]),
            height=height[i],
            radius=radius[i],
        )
        phase = Phase(shape=elem.generate())
        phases.append(phase)

    return [periodic(phase=phase_elem, rve=rve) for phase_elem in phases]


@pytest.fixture()
def box_homogeneous_unit(rve_unit: Rve) -> tuple[cq.Shape, list[Phase], Rve]:
    """Return a homogeneous unit box."""
    shape = Box(
        center=tuple(rve_unit.center),
        orientation=(0.0, 0.0, 0.0),
        dim=(rve_unit.dim[0], rve_unit.dim[1], rve_unit.dim[2]),
    ).generate()
    listcqphases = [Phase(shape=shape)]
    return (shape, listcqphases, rve_unit)


@pytest.fixture()
def box_homogeneous_double(
    rve_double: Rve,
) -> tuple[cq.Shape, list[Phase], Rve]:
    """Return a homogeneous double box."""
    shape = Box(
        center=tuple(rve_double.center),
        orientation=(0.0, 0.0, 0.0),
        dim=tuple(rve_double.dim),
    ).generate()
    listcqphases = [Phase(shape=shape)]
    return (shape, listcqphases, rve_double)


@pytest.fixture()
def box_homogeneous_double_centered(
    rve_double_centered: Rve,
) -> (cq.Shape, list[Phase], Rve):
    """Return a homogeneous double centered box."""
    shape = Box(
        center=tuple(rve_double_centered.center),
        orientation=(0.0, 0.0, 0.0),
        dim=tuple(rve_double_centered.dim),
    ).generate()
    listcqphases = [Phase(shape=shape)]
    return (shape, listcqphases, rve_double_centered)


@pytest.fixture()
def octet_truss_homogeneous_unit(
    rve_unit: Rve,
) -> tuple[cq.Shape, list[Phase], Rve]:
    """Return a homogeneous unit octet truss."""
    list_periodic_phases = _generate_cqcompound_octettruss(rve_unit)
    merged = fuseShapes(
        [phase.shape for phase in list_periodic_phases],
        retain_edges=False,
    )
    listcqphases = [Phase(shape=merged)]
    return (merged, listcqphases, rve_unit)


@pytest.fixture()
def octet_truss_homogeneous_double_centered(
    rve_double_centered: Rve,
) -> tuple[cq.Shape, list[Phase], Rve]:
    """Return a homogeneous double centered octet truss."""
    list_periodic_phases = _generate_cqcompound_octettruss(rve_double_centered)
    merged = fuseShapes(
        [phase.shape for phase in list_periodic_phases],
        retain_edges=False,
    )
    listcqphases = [Phase(shape=merged)]
    return (merged, listcqphases, rve_double_centered)


@pytest.fixture()
def octet_truss_heterogeneous(
    rve_unit: Rve,
) -> tuple[cq.Compound, list[Phase], Rve]:
    """Return a heterogeneous octet truss."""
    list_periodic_phases = _generate_cqcompound_octettruss(rve_unit)
    listcqphases = cutPhases(phaseList=list_periodic_phases, reverseOrder=False)
    return (
        cq.Compound.makeCompound([phase.shape for phase in listcqphases]),
        listcqphases,
        rve_unit,
    )


@pytest.mark.filterwarnings("ignore:Object intersecting")
@pytest.mark.parametrize(
    "shape",
    [
        "box_homogeneous_unit",
        "box_homogeneous_double",
        "box_homogeneous_double_centered",
        "octet_truss_homogeneous_unit",
        "octet_truss_homogeneous_double_centered",
        "octet_truss_heterogeneous",
    ],
)
def test_octettruss_mesh_must_be_periodic(
    shape: tuple[cq.Compound | cq.Shape, list[Phase], Rve],
    request: pytest.FixtureRequest,
    tmp_output_compound_filename: str,
    tmp_output_vtk_filename: str,
) -> None:
    """Test that the octet truss mesh is periodic."""
    cqoctet, listcqphases, rve = request.getfixturevalue(shape)

    cq.exporters.export(cqoctet, tmp_output_compound_filename)
    meshPeriodic(
        mesh_file=tmp_output_compound_filename,
        rve=rve,
        listPhases=listcqphases,
        size=0.03,
        order=1,
        output_file=tmp_output_vtk_filename,
    )

    # Act
    pvmesh = pv.read(tmp_output_vtk_filename)
    crd = pvmesh.points

    assert is_periodic(crd)
