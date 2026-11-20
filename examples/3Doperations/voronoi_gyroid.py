"""
Voronoi gyroid
==============

Intersect each cell of a Voronoi tessellation with a gyroid sheet, then mesh
the resulting phases with gmsh.

Requires the CAD backend: ``pip install microgen[cad]``.
"""

from microgen import Neper, Phase, Tpms, mesh
from microgen.cad import make_compound
from microgen.shape import surface_functions

tess_file = "voronoi.tess"
polyhedra = Neper.voronoi_from_tess_file(tess_file)

gyroid = Tpms(
    center=(0.5, 0.5, 0.5),
    surface_function=surface_functions.gyroid,
    offset=0.2,
)
gyroid = gyroid.generate_cad(type_part="sheet").translate((0.5, 0.5, 0.5))

phases = []
for polyhedron in polyhedra:
    shape = polyhedron.generate_cad()
    phases.append(Phase.from_cad(shape.intersect(gyroid)))

compound = make_compound([phase.cad for phase in phases])
step_file = "voronoi_gyroid.step"
compound.export_step(step_file)

vtk_file = "voronoi_gyroid.vtk"
mesh(
    mesh_file=step_file,
    list_phases=phases,
    size=0.05,
    order=1,
    output_file=vtk_file,
)

# %%
# .. image:: /_static/examples/Gyroid-voro.png
#
# sphinx_gallery_thumbnail_path = "_static/examples/Gyroid-voro.png"
