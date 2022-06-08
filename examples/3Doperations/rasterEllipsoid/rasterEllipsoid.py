import cadquery as cq
from microgen import Rve, Ellipsoid, rasterShapeList, mesh, Phase

rve = Rve(dim_x=1, dim_y=1, dim_z=1)
listPhases = []

elem = Ellipsoid(center=(0, 0.5, 0.5),
                 a_x=0.15, a_y=0.31, a_z=0.4)
elli = elem.generate()

raster = rasterShapeList(cqShapeList=[elli], rve=rve, grid=[5, 5, 5])

phases = [Phase(solids=solids) for solids in raster[1]]
flat_list = [solid for phase in phases for solid in phase.getSolids()]

compound = cq.Compound.makeCompound(flat_list)
cq.exporters.export(compound, 'compound.step')

mesh(mesh_file='compound.step', listPhases=phases, size=0.03, order=1, output_file='rasterEllipsoid.vtk')
