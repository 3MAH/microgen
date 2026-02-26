"""
This script generates a gyroid-based TPMS (Triply Periodic Minimal Surface) structure
using *cadquery*, exports the geometry to a STEP file, performs periodic
meshing, and optionally remeshes the result for finite element method (FEM) simulations.

Steps performed:
1. Generate a TPMS geometry using the gyroid surface function.
2. Wrap the geometry into a microgen Phase object.
3. Export the geometry as a STEP file.
4. Import the STEP file and create a mesh with periodic constraints and export as VTK.
5. Optionally remesh the mesh while preserving periodicity.

Generated files:
- gyroid.step: CAD representation of the TPMS geometry.
- gyroid.vtk: Initial meshed structure.
- remeshed_gyroid_mesh.vtk (optional): Quality-improved periodic mesh for FEM.

Dependencies:
- microgen
- cadquery
- pyvista
"""

from pathlib import Path
from microgen import Tpms, Phase, meshPeriodic, Rve
from microgen.remesh import remesh_keeping_periodicity_for_fem
from microgen.shape.surface_functions import gyroid
import cadquery as cq
import pyvista as pv

# 1. Generate a TPMS geometry using the gyroid surface function.
geometry = Tpms(
    surface_function=gyroid,
    density=0.30,
    resolution=30,
)

# 2. Wrap the geometry into a microgen Phase object.
phases = []
phases.append(Phase(shape=geometry.generate()))
rve = Rve(dim=1)

# 3. Export the geometry as a STEP file.
step_file = str(Path(__file__).parent / "gyroid.step")
cq.exporters.export(phases[0].shape, step_file)

# 4. Import the STEP file and create a mesh with periodic constraints and export as VTK.

vtk_file = str(Path(__file__).parent / "gyroid.vtk")

meshPeriodic(
    mesh_file=step_file,
    rve=rve,
    listPhases=phases,
    order=1,
    size=0.03,
    output_file=vtk_file,
)

# 5. Optionally remesh the mesh while preserving periodicity.

initial_gyroid = pv.UnstructuredGrid(vtk_file)
max_element_edge_length = 0.02
remeshed_gyroid = remesh_keeping_periodicity_for_fem(
    initial_gyroid,
    hmax=max_element_edge_length,
)
remeshed_vtk_file = str(Path(__file__).parent / "remeshed_gyroid_mesh.vtk")
# remeshed_gyroid.save(remeshed_vtk_file)
