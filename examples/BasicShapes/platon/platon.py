import microgen
import cadquery as cq


filenames = [
    "tetrahedron.obj",
    "cube.obj",
    "octahedron.obj",
    "dodecahedron.obj",
    "icosahedron.obj",
]
platon_solids = []

i = 0
for filename in filenames:
    dic = microgen.shape.polyhedron.read_obj(filename)
    solid = microgen.shape.polyhedron.Polyhedron(dic=dic)
    shape = solid.generate()
    shape = shape.translate(cq.Vector(-6 + 3 * i, 0, 0))
    platon_solids.append(shape)
    i += 1

cq.exporters.export(cq.Compound.makeCompound(platon_solids), "platon.stl")

