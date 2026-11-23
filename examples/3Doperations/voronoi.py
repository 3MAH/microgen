"""
Voronoi polycrystal
===================

Import the cells of a Neper tessellation as polyhedra with
:meth:`microgen.Neper.voronoi_from_tess_file`, build a CAD compound and mesh
it with gmsh.

Requires the CAD backend: ``pip install microgen[cad]``.
"""

import pyvista as pv

from microgen import Neper, Phase, mesh
from microgen.cad import make_compound

tess_file = "voronoi.tess"
polyhedra = Neper.voronoi_from_tess_file(tess_file)

shapes = [poly.generate_cad() for poly in polyhedra]

compound = make_compound(shapes)
step_file = "voronoi.step"
compound.export_step(step_file)

phases = [Phase.from_cad(shape) for shape in shapes]

vtk_file = "voronoi.vtk"
mesh(
    mesh_file=step_file,
    list_phases=phases,
    size=0.05,
    order=1,
    output_file=vtk_file,
)

# %%
# One color per grain.
grid = pv.read(vtk_file)
grains = grid.cell_data.keys()[0]
grid.plot(
    scalars=grains,
    cmap="tab20",
    smooth_shading=True,
    split_sharp_edges=True,
    show_scalar_bar=False,
)

# %%
# The tetrahedral mesh: each grain is a separate physical volume, with
# matching nodes on the grain boundaries.
grid.plot(show_edges=True, scalars=grains, cmap="tab20", show_scalar_bar=False)
