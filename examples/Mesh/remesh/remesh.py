import os
import warnings
from distutils.spawn import find_executable

import pyvista as pv

from microgen import Tpms
from microgen.remesh import remesh_keeping_periodicity_for_fem
from microgen.shape.surface_functions import gyroid

USE_MMG = find_executable("mmg3d_O3") is not None
if not USE_MMG:
    warnings.warn("MMG not found")

if USE_MMG:
    if "data" not in os.listdir("."):
        os.mkdir("data")

    print("generate gyroid", flush=True)
    initial_gyroid = pv.UnstructuredGrid(
        Tpms(surface_function=gyroid, offset=1.0, resolution=50).generateVtk(
            type_part="sheet"
        )
    )
    print("save gyroid", flush=True)
    initial_gyroid.save("data/initial_gyroid_mesh.vtk")

    print("remesh gyroid", flush=True)
    max_element_edge_length = 0.02
    remeshed_gyroid = remesh_keeping_periodicity_for_fem(
        initial_gyroid, hmax=max_element_edge_length
    )
    print("save remeshed gyroid", flush=True)
    remeshed_gyroid.save("data/remeshed_gyroid_mesh.vtk")
