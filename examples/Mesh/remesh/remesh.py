import os

import pyvista as pv

from microgen import BoxMesh, Tpms
from microgen.remesh import remesh_keeping_periodicity_for_fem
from microgen.shape.surface_functions import gyroid

if "data" not in os.listdir("."):
    os.mkdir("data")

gyroid_vtk = pv.UnstructuredGrid(
    Tpms(surface_function=gyroid, offset=1.0, resolution=50).generateVtk(
        type_part="sheet"
    )
)
gyroid_vtk.save("data/initial_gyroid_mesh.vtk")
gyroid_mesh = BoxMesh.from_pyvista(gyroid_vtk)

remeshed_gyroid_mesh = remesh_keeping_periodicity_for_fem(gyroid_mesh, hmax=0.02)
remeshed_gyroid_mesh.save("data/remeshed_gyroid_mesh.vtk")
