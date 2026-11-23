"""
Preset strut lattices
=====================

The twelve preset strut-based lattices of :mod:`microgen.shape.strut_lattice`,
with a strut radius of 0.1.

Each lattice is turned into one implicit shape with
:meth:`~microgen.AbstractLattice.generate_implicit` (union of the struts and
joints, clipped by the unit cell) and meshed by marching cubes, without the
CAD backend.
"""

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
    BodyCenteredCubic(strut_radius=0.1),
    Cubic(strut_radius=0.1),
    Cuboctahedron(strut_radius=0.1),
    Diamond(strut_radius=0.1),
    FaceCenteredCubic(strut_radius=0.1),
    Octahedron(strut_radius=0.1),
    OctetTruss(strut_radius=0.1),
    RhombicCuboctahedron(strut_radius=0.1),
    RhombicDodecahedron(strut_radius=0.1),
    TruncatedCube(strut_radius=0.1),
    TruncatedCuboctahedron(strut_radius=0.1),
    TruncatedOctahedron(strut_radius=0.1),
]

# %%
# One subplot per lattice, on a 3 x 4 grid.
bounds = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
N_COL = 4
pl = pv.Plotter(shape=(3, N_COL), window_size=(1200, 900), border=False)
for i, lattice in enumerate(preset_lattice_list):
    mesh = lattice.generate_implicit().generate_surface_mesh(
        bounds=bounds, resolution=80
    )
    pl.subplot(i // N_COL, i % N_COL)
    pl.add_text(type(lattice).__name__, font_size=8)
    pl.add_mesh(mesh, color="#ef7d2d", smooth_shading=True)
pl.link_views()
pl.show()
