from pathlib import Path

import cadquery as cq

import microgen

filenames = [
    "tetrahedron.obj",
    "cube.obj",
    "octahedron.obj",
    "dodecahedron.obj",
    "icosahedron.obj",
]
filenames = [str(Path(__file__).parent / file) for file in filenames]
platon_solids = []

i = 0
for filename in filenames:
    dic = microgen.shape.polyhedron.read_obj(filename)
    solid = microgen.shape.polyhedron.Polyhedron(dic=dic)
    shape = solid.generate()
    shape = shape.translate(cq.Vector(-6 + 3 * i, 0, 0))
    platon_solids.append(shape)
    i += 1

stl_file = str(Path(__file__).parent / "platon.stl")
cq.exporters.export(cq.Compound.makeCompound(platon_solids), stl_file)
