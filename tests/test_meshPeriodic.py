from pathlib import Path
from typing import Union

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


@pytest.fixture(scope="function")
def rve_unit() -> Rve:
    return Rve(dim=1, center=(0.5, 0.5, 0.5))


@pytest.fixture(scope="function")
def rve_double() -> Rve:
    return Rve(dim=2, center=(1.0, 1.0, 1.0))


@pytest.fixture(scope="function")
def rve_double_centered() -> Rve:
    return Rve(dim=2, center=(0.0, 0.0, 0.0))


@pytest.fixture(scope="function")
def tmp_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    tmp_dir_name = "test_tmp_dir"
    return tmp_path_factory.mktemp(tmp_dir_name)


@pytest.fixture(scope="function")
def tmp_output_compound_filename(tmp_dir: Path) -> str:
    return (tmp_dir / "compound.step").as_posix()


@pytest.fixture(scope="function")
def tmp_output_vtk_filename(tmp_dir: Path) -> str:
    return (tmp_dir / "octettruss.vtk").as_posix()


@pytest.mark.filterwarnings("ignore:Object intersecting")
def _generate_cqcompound_octettruss(rve: Rve):
    # data to generate octet truss cylinders
    xc = np.array([0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0])
    yc = np.array([0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5])
    zc = np.array([0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5])
    psi = np.array(
        [45.0, -45.0, 0.0, 0.0, 90.0, 90.0, 0.0, 0.0, 90.0, 90.0, 45.0, -45.0]
    )
    theta = np.array(
        [0.0, 0.0, -90.0, -90.0, -90.0, -90.0, 90.0, 90.0, 90.0, 90.0, 0.0, 0.0]
    )
    phi = np.array(
        [0.0, 0.0, 45.0, -45.0, 45.0, -45.0, -45.0, 45.0, -45.0, 45.0, 0.0, 0.0]
    )
    height = np.full_like(xc, 2.0)
    radius = np.full_like(xc, 0.05)

    listPhases: list[Phase] = []
    listPeriodicPhases: list[Phase] = []
    n = len(xc)

    for i in range(0, n):
        elem = Cylinder(
            center=(xc[i], yc[i], zc[i]),
            orientation=(psi[i], theta[i], phi[i]),
            height=height[i],
            radius=radius[i],
        )
        phase = Phase(shape=elem.generate())
        listPhases.append(phase)

    for phase_elem in listPhases:
        periodicPhase = periodic(phase=phase_elem, rve=rve)
        listPeriodicPhases.append(periodicPhase)

    return listPeriodicPhases


@pytest.fixture(scope="function")
def box_homogeneous_unit(rve_unit: Rve) -> (cq.Shape, list[Phase]):
    shape = Box(
        center=rve_unit.center,
        orientation=(0.0, 0.0, 0.0),
        dim_x=rve_unit.dim[0],
        dim_y=rve_unit.dim[1],
        dim_z=rve_unit.dim[2],
    ).generate()
    listcqphases = [Phase(shape=shape)]
    return (shape, listcqphases, rve_unit)


@pytest.fixture(scope="function")
def box_homogeneous_double(rve_double: Rve) -> (cq.Shape, list[Phase]):
    shape = Box(
        center=rve_double.center,
        orientation=(0.0, 0.0, 0.0),
        dim_x=rve_double.dim[0],
        dim_y=rve_double.dim[1],
        dim_z=rve_double.dim[2],
    ).generate()
    listcqphases = [Phase(shape=shape)]
    return (shape, listcqphases, rve_double)


@pytest.fixture(scope="function")
def box_homogeneous_double_centered(
    rve_double_centered: Rve,
) -> (cq.Shape, list[Phase]):
    shape = Box(
        center=rve_double_centered.center,
        orientation=(0.0, 0.0, 0.0),
        dim_x=rve_double_centered.dim[0],
        dim_y=rve_double_centered.dim[1],
        dim_z=rve_double_centered.dim[2],
    ).generate()
    listcqphases = [Phase(shape=shape)]
    return (shape, listcqphases, rve_double_centered)


@pytest.fixture(scope="function")
def octet_truss_homogeneous_unit(rve_unit: Rve) -> (cq.Shape, list[Phase]):
    listPeriodicPhases = _generate_cqcompound_octettruss(rve_unit)
    merged = fuseShapes(
        [phase.shape for phase in listPeriodicPhases], retain_edges=False
    )
    listcqphases = [Phase(shape=merged)]
    return (merged, listcqphases, rve_unit)


@pytest.fixture(scope="function")
def octet_truss_homogeneous_double_centered(
    rve_double_centered: Rve,
) -> (cq.Shape, list[Phase]):
    listPeriodicPhases = _generate_cqcompound_octettruss(rve_double_centered)
    merged = fuseShapes(
        [phase.shape for phase in listPeriodicPhases], retain_edges=False
    )
    listcqphases = [Phase(shape=merged)]
    return (merged, listcqphases, rve_double_centered)


@pytest.fixture(scope="function")
def octet_truss_heterogeneous(rve_unit: Rve) -> (cq.Compound, list[Phase]):
    listPeriodicPhases = _generate_cqcompound_octettruss(rve_unit)
    listcqphases = cutPhases(phaseList=listPeriodicPhases, reverseOrder=False)
    return (
        cq.Compound.makeCompound([phase.shape for phase in listcqphases]),
        listcqphases,
        rve_unit,
    )


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
    shape: Union[cq.Compound, cq.Shape],
    request,
    tmp_output_compound_filename: str,
    tmp_output_vtk_filename: str,
):
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
