"""
Periodic gyroid mesh from a STEP file
=====================================

Generate a gyroid, export it to a STEP file, mesh it periodically with gmsh,
and remesh the result with mmg for finite element simulations.

Steps performed:

1. Generate a TPMS geometry using the gyroid surface function.
2. Wrap the geometry into a microgen Phase object.
3. Export the geometry as a STEP file.
4. Import the STEP file and create a mesh with periodic constraints and export as VTK.
5. Remesh the mesh while preserving periodicity.

Generated files:

- periodic_gyroid.step: CAD representation of the TPMS geometry.
- periodic_gyroid.vtk: Initial meshed structure.
- periodic_gyroid_remeshed.vtk: Quality-improved periodic mesh for FEM.

Requires the CAD backend (``pip install microgen[cad]``) and the ``mmg``
executables (``conda install -c conda-forge mmg``).
"""

import pyvista as pv

from microgen import Phase, Rve, Tpms, mesh_periodic
from microgen.remesh import remesh_keeping_boundaries_for_fem
from microgen.shape.surface_functions import gyroid

# 1. Generate a TPMS geometry using the gyroid surface function.
geometry = Tpms(
    surface_function=gyroid,
    density=0.30,
    resolution=30,
)

# 2. Wrap the geometry into a microgen Phase object.
phases = []
phases.append(Phase.from_cad(geometry.generate_cad()))
rve = Rve(dim=1)

# 3. Export the geometry as a STEP file.
step_file = "periodic_gyroid.step"
phases[0].cad.export_step(step_file)

# 4. Import the STEP file and create a mesh with periodic constraints and export as VTK.

vtk_file = "periodic_gyroid.vtk"

mesh_periodic(
    mesh_file=step_file,
    rve=rve,
    list_phases=phases,
    order=1,
    size=0.03,
    output_file=vtk_file,
)

# 5. Remesh the mesh while preserving periodicity.

initial_gyroid = pv.UnstructuredGrid(vtk_file)
max_element_edge_length = 0.02
remeshed_gyroid = remesh_keeping_boundaries_for_fem(
    initial_gyroid,
    periodic=True,
    hmax=max_element_edge_length,
)
remeshed_vtk_file = "periodic_gyroid_remeshed.vtk"
remeshed_gyroid.save(remeshed_vtk_file)

# %%
# .. image:: /_static/examples/remeshed_mesh.png
#
# sphinx_gallery_thumbnail_path = "_static/examples/remeshed_mesh.png"
