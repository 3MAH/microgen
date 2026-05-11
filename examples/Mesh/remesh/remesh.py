"""Remesh a gyroid surface keeping periodicity for FEM simulations."""

from pathlib import Path

from microgen import Tpms
from microgen.remesh import remesh_keeping_boundaries_for_fem
from microgen.shape.surface_functions import gyroid

data_dir = Path(__file__).parent / "data"
Path.mkdir(data_dir, exist_ok=True)

tpms = Tpms(surface_function=gyroid, offset=1.0, resolution=50)
initial_gyroid = tpms.grid_sheet
initial_gyroid.save(data_dir / "initial_gyroid_mesh.vtk")

max_element_edge_length = 0.02
remeshed_gyroid = remesh_keeping_boundaries_for_fem(
    initial_gyroid,
    hmax=max_element_edge_length,
)
remeshed_gyroid.save(data_dir / "remeshed_gyroid_mesh.vtk")
