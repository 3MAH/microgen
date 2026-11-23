"""
Platonic solids
===============

Read the five Platonic solids from Wavefront ``.obj`` files with
:func:`microgen.shape.polyhedron.read_obj` and build them as
:class:`microgen.Polyhedron` shapes.

Requires the CAD backend: ``pip install microgen[cad]``.
"""

import pyvista as pv

import microgen
from microgen.cad import make_compound

filenames = [
    "tetrahedron.obj",
    "cube.obj",
    "octahedron.obj",
    "dodecahedron.obj",
    "icosahedron.obj",
]

platon_solids = []
for i, filename in enumerate(filenames):
    dic = microgen.shape.polyhedron.read_obj(filename)
    solid = microgen.shape.polyhedron.Polyhedron(dic=dic)
    shape = solid.generate_cad()
    shape = shape.translate((-6 + 3 * i, 0, 0))
    platon_solids.append(shape)

stl_file = "platonic_solids.stl"
make_compound(platon_solids).export_stl(stl_file)

# %%
pv.read(stl_file).plot(show_edges=True, color="#ef7d2d")
