from pathlib import Path
from microgen import Tpms, Phase, meshPeriodic, Rve
from microgen.remesh import remesh_keeping_periodicity_for_fem
from microgen.shape.surface_functions import gyroid
import cadquery as cq
import pyvista as pv

geometry = Tpms(
    surface_function=gyroid,
    density=0.30,
    resolution=30,
)

listPhases = []
listPhases.append(Phase(shape=geometry.generate()))
rve = Rve(dim=1)

step_file = str(Path(__file__).parent / "gyroid.step")
cq.exporters.export(listPhases[0].shape, step_file)

vtk_file = str(Path(__file__).parent / "gyroid.vtk")

meshPeriodic(
    mesh_file=step_file,
    rve=rve,
    listPhases=listPhases,
    order=1,
    size=0.03,
    output_file=vtk_file,
)

initial_gyroid = pv.UnstructuredGrid(vtk_file)

max_element_edge_length = 0.02
remeshed_gyroid = remesh_keeping_periodicity_for_fem(
    initial_gyroid,
    hmax=max_element_edge_length,
)
remeshed_vtk_file = str(Path(__file__).parent / "remeshed_gyroid_mesh.vtk")
# remeshed_gyroid.save(remeshed_vtk_file)
