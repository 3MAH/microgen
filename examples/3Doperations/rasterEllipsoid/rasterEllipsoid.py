import cadquery as cq
from microgen import Rve, Ellipsoid, rasterShapeList, mesh, Phase

rve = Rve(dim_x=1, dim_y=1, dim_z=1)
listPhases = []

elem = Ellipsoid(center=(0, 0.5, 0.5),
                 a_x=0.15, a_y=0.31, a_z=0.4)
elli = elem.generate()

cq.exporters.export(elli, 'ellipsoid.step')
raster = rasterShapeList(cqShapeList=[elli], rve=rve, grid=[5, 5, 5])
phase = Phase(solids=raster[1])

compound = cq.Compound.makeCompound(phase.getFlatSolidList())
cq.exporters.export(compound, 'compound.step')

mesh(mesh_file='compound.step', listPhases=[phase], size=0.03, order=1, output_file='Mesh.msh')
