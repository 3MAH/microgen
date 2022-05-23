import cadquery as cq
from microgen import Rve, BasicGeometry, rasterShapeList, mesh, gyroid

# Size of the mesh
size_mesh = 0.03
a = 1.0
b = 1.0
c = 1.0

periodicity = 0
revel = Rve(a, b, c, size_mesh)

elem = BasicGeometry(number=0, shape='tpms',
                     xc=0.5, yc=0.5, zc=0.5,
                     psi=0., theta=0., phi=0.,
                     param_geom={"surface_function": gyroid,
                                 "type_part": "skeletal",
                                 "thickness": 0.1},
                     path_data='data')
elem.geometry.createSurfaces(rve=revel,
                             sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                             path_data='data')
skeletal = elem.generate(rve=revel)

cq.exporters.export(skeletal, 'skeletal.step')

raster = rasterShapeList(cqShapeList=[skeletal], rve=revel, grid=[5, 5, 5])

compound = cq.Compound.makeCompound(raster[0])
cq.exporters.export(compound, 'compound.step')

mesh(mesh_file='compound.step', listPhases=raster[1], size=size_mesh, order=1, output_file='Mesh.msh')
