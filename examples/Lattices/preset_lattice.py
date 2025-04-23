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
    BodyCenteredCubic(),
    Cubic(),
    Cuboctahedron(),
    Diamond(),
    FaceCenteredCubic(),
    Octahedron(),
    OctetTruss(),
    RhombicCuboctahedron(),
    RhombicDodecahedron(),
    TruncatedCube(),
    TruncatedCuboctahedron(),
    TruncatedOctahedron(),
]

meshes = []

i = 0
n_col = 4
n_row = np.ceil(len(preset_lattice_list) / n_col)
for i, lattice in enumerate(preset_lattice_list):
    i_x = i % n_col
    i_y = i // n_col
    mesh = lattice.generate_vtk()
    mesh.translate(
        [1.2 * (i_x - 0.5 * (n_col - 1)), -1.2 * (i_y - 0.5 * (n_row - 1)), 0],
        inplace=True,
    )
    meshes.append(mesh)

vtm_file = Path(__file__).parent / "lattices.vtm"
pv.MultiBlock(meshes).save(vtm_file)
