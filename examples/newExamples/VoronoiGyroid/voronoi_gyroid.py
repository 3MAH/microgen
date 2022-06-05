import cadquery as cq
from microgen import parseNeper, Rve, BasicGeometry, mesh, cutPhaseByShapeList
from microgen.operations import cutPhasesByShape
import microgen

# We import the Polyhedra from Nepoer tesselation file
listPolyhedra, seed, vertices, edges, faces, polys = parseNeper('test1')

# We fix the dimensions of the RVE
a = 1.0
b = 1.0
c = 1.0

rve = Rve(a, b, c)

gyroid = BasicGeometry(number=0, shape='tpms',
                       xc=0.5, yc=0.5, zc=0.5,
                       psi=0., theta=0., phi=0.,
                       param_geom={"surface_function": microgen.shape.tpms.gyroid,
                                   "type_part": 'sheet',
                                   "thickness": 0.2},
                       path_data='data')
gyroid.geometry.createSurfaces(rve=rve,
                               sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                               path_data='data')
gyroid = gyroid.generate(rve=rve).translate((0.5, 0.5, 0.5))

phases = []
i = 0
for polyhedron in listPolyhedra:
    elem = BasicGeometry(number=0, shape='polyhedron',
                         xc=polyhedron["original"][0],
                         yc=polyhedron["original"][1],
                         zc=polyhedron["original"][2],
                         psi=0., theta=0., phi=0.,
                         param_geom={"dic": polyhedron}, path_data='.')
    elem = elem.generate()
    compound = elem.intersect(gyroid)
    cq.exporters.export(compound, 'compound{}.stl'.format(i))
    i = i + 1
