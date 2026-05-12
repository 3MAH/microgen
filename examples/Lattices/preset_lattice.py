from pathlib import Path

import numpy as np
import pyvista as pv

from microgen import (
    BodyCenteredCubic,
    Cubic,
    Cuboctahedron,
    Diamond,
    FaceCenteredCubic,
    Octahedron,
    OctetTruss,
    RhombicCuboctahedron,
    RhombicDodecahedron,
    TruncatedCube,
    TruncatedCuboctahedron,
    TruncatedOctahedron,
)

preset_lattice_list = [
    BodyCenteredCubic(0.1),
    Cubic(0.1),
    Cuboctahedron(0.1),
    Diamond(0.1),
    FaceCenteredCubic(0.1),
    Octahedron(0.1),
    OctetTruss(0.1),
    RhombicCuboctahedron(0.1),
    RhombicDodecahedron(0.1),
    TruncatedCube(0.1),
    TruncatedCuboctahedron(0.1),
    TruncatedOctahedron(0.1),
]

meshes = pv.PolyData()

N_COL = 4
n_row = np.ceil(len(preset_lattice_list) / N_COL)
for i, lattice in enumerate(preset_lattice_list):
    i_x = i % N_COL
    i_y = i // N_COL
    mesh = lattice.generate_surface_mesh()
    mesh.translate(
        [1.2 * (i_x - 0.5 * (N_COL - 1)), -1.2 * (i_y - 0.5 * (n_row - 1)), 0],
        inplace=True,
    )
    meshes.append_polydata(mesh, inplace=True)

stl_file = Path(__file__).parent / "lattices.stl"
meshes.save(stl_file)
