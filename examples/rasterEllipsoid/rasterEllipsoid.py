import cadquery as cq
from microgen import Rve, BasicGeometry, rasterShapeList, mesh

# Size of the mesh
size_mesh = 0.03
a = 1.0
b = 1.0
c = 1.0

periodicity = 0
revel = Rve(a, b, c, size_mesh)
listPhases = []

elem = BasicGeometry(number=101, shape='ellipsoid',
                     xc=0., yc=0.5, zc=0.5,
                     psi=0, theta=0, phi=0,
                     param_geom={"a_x": 0.15,
                                 "a_y": 0.31,
                                 "a_z": 0.4},
                     path_data='')
elli = elem.generate(rve=revel)

cq.exporters.export(elli, 'ellipsoid.step')
raster = rasterShapeList(cqShapeList=[elli], rve=revel, grid=[5, 5, 5])

compound = cq.Compound.makeCompound(raster[0])
cq.exporters.export(compound, 'compound.step')

mesh(mesh_file='compound.step', listPhases=raster[1], size=size_mesh, order=1, output_file='Mesh.msh')
