import os
from pathlib import Path

import pyvista as pv

from microgen import BoxMesh, Tpms
from microgen.remesh import remesh_keeping_periodicity_for_fem
from microgen.shape.surface_functions import gyroid

data_dir = Path(__file__).parent / "data"
os.makedirs(data_dir, exist_ok=True)

print("generate gyroid", flush=True)
gyroid_vtk = pv.UnstructuredGrid(
    Tpms(surface_function=gyroid, offset=1.0, resolution=50).generateVtk(
        type_part="sheet"
    )
)
print("save gyroid", flush=True)
gyroid_vtk.save(data_dir / "initial_gyroid_mesh.vtk")
print("convert gyroid", flush=True)
gyroid_mesh = BoxMesh.from_pyvista(gyroid_vtk)

print("remesh gyroid", flush=True)
remeshed_gyroid_mesh = remesh_keeping_periodicity_for_fem(gyroid_mesh, hmax=0.02)
print("save remeshed gyroid", flush=True)
remeshed_gyroid_mesh.save(data_dir / "remeshed_gyroid_mesh.vtk")
