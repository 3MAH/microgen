import os
from pathlib import Path

import pyvista as pv

from microgen import Tpms
from microgen.remesh import remesh_keeping_periodicity_for_fem
from microgen.shape.surface_functions import gyroid

data_dir = Path(__file__).parent / "data"
os.makedirs(data_dir, exist_ok=True)

print("generate gyroid", flush=True)
initial_gyroid = pv.UnstructuredGrid(
    Tpms(surface_function=gyroid, offset=1.0, resolution=50).generate_vtk(
        type_part="sheet",
    ),
)
print("save gyroid", flush=True)
initial_gyroid.save(data_dir / "initial_gyroid_mesh.vtk")

print("remesh gyroid", flush=True)
max_element_edge_length = 0.02
remeshed_gyroid = remesh_keeping_periodicity_for_fem(
    initial_gyroid,
    hmax=max_element_edge_length,
)
print("save remeshed gyroid", flush=True)
remeshed_gyroid.save(data_dir / "remeshed_gyroid_mesh.vtk")
