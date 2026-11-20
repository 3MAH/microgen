"""
Periodic remeshing
==================

Remesh a gyroid volume mesh with mmg while keeping its boundaries, so that
the mesh stays periodic for finite element simulations.

Requires the ``mmg`` executables: ``conda install -c conda-forge mmg``.
"""

from pathlib import Path

from microgen import Tpms
from microgen.remesh import remesh_keeping_boundaries_for_fem
from microgen.shape.surface_functions import gyroid

data_dir = Path("remesh_data")
data_dir.mkdir(exist_ok=True)

tpms = Tpms(surface_function=gyroid, offset=1.0, resolution=50)
initial_gyroid = tpms.grid_sheet
initial_gyroid.save(data_dir / "initial_gyroid_mesh.vtk")

max_element_edge_length = 0.02
remeshed_gyroid = remesh_keeping_boundaries_for_fem(
    initial_gyroid,
    hmax=max_element_edge_length,
)
remeshed_gyroid.save(data_dir / "remeshed_gyroid_mesh.vtk")

# %%
# Initial mesh:
#
# .. image:: /_static/examples/initial_mesh.png
#
# Remeshed mesh:
#
# .. image:: /_static/examples/remeshed_mesh.png
#
# sphinx_gallery_thumbnail_path = "_static/examples/remeshed_mesh.png"
