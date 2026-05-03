from pathlib import Path

import microgen
from microgen.cad import make_compound

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
    shape = shape.translate((-6 + 3 * i, 0, 0))
    platon_solids.append(shape)
    i += 1

stl_file = str(Path(__file__).parent / "platon.stl")
make_compound(platon_solids).export_stl(stl_file)
