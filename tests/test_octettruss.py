from microgen import Rve, Phase, Cylinder, periodic, fuseShapes, cutPhases, meshPeriodic, is_periodic
import cadquery as cq
import pyvista as pv
import numpy as np
from typing import Union
import pytest
from pathlib import Path

@pytest.fixture(scope='function')
def rve() -> Rve:
    return Rve(dim_x=1, dim_y=1, dim_z=1, center=(0.5, 0.5, 0.5))

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
def _generate_cqcompound_octettruss(
    rve: Rve
):

    #data to generate octet truss cylinders
    xc = np.array([0.5,0.5,0.5,0.5,0.0,0.0,0.5,0.5,0.5,0.5,0.0,0.0])
    yc = np.array([0.5,0.5,0.0,0.0,0.5,0.5,0.5,0.5,0.0,0.0,0.5,0.5])
    zc = np.array([0.0,0.0,0.5,0.5,0.5,0.5,0.0,0.0,0.5,0.5,0.5,0.5])
    psi = np.array([45.,-45.,0.,0.,90.,90.,0.,0.,90.,90.,45.,-45.])
    theta = np.array([0.,0.,-90.,-90.,-90.,-90.,90.,90.,90.,90.,0.,0.])
    phi = np.array([0.,0.,45.,-45.,45.,-45.,-45.,45.,-45.,45.,0.,0.])
    height = np.full_like(xc, 2.0)
    radius = np.full_like(xc, 0.05)

    listPhases : list[Phase] = []
    listPeriodicPhases : list[Phase] = []
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
def octet_truss_homogeneous(rve : Rve) -> (cq.Shape, list[Phase]):

    listPeriodicPhases = _generate_cqcompound_octettruss(rve)
    merged = fuseShapes([phase.shape for phase in listPeriodicPhases], retain_edges=False)
    listcqphases = [Phase(shape=merged)]    
    return (merged, listcqphases)

@pytest.fixture(scope="function")
def octet_truss_heterogeneous(rve : Rve) -> (cq.Compound, list[Phase]):

    listPeriodicPhases = _generate_cqcompound_octettruss(rve)    
    listcqphases = cutPhases(phaseList=listPeriodicPhases, reverseOrder=False)
    return (cq.Compound.makeCompound([phase.shape for phase in listcqphases]), listcqphases)

@pytest.mark.parametrize(
    "shape",
    [
        "octet_truss_homogeneous",
        "octet_truss_heterogeneous",
    ]
)
def test_octettruss_mesh_must_be_periodic(
    shape: Union[cq.Compound, cq.Shape],
    request, 
    rve: Rve,
    tmp_output_compound_filename: str,
    tmp_output_vtk_filename: str
):

    cqoctet, listcqphases = request.getfixturevalue(shape)

    cq.exporters.export(cqoctet, tmp_output_compound_filename)
    meshPeriodic(
        mesh_file=tmp_output_compound_filename,
        rve=rve,
        listPhases=listcqphases,
        size=0.03,
        order=1,
        output_file=tmp_output_vtk_filename,
    )

    #Act
    pvmesh = pv.read(tmp_output_vtk_filename)  
    crd = pvmesh.points

    assert is_periodic(crd)    
