import cadquery as cq

from microgen import Ellipsoid, Phase, Rve, mesh, rasterPhase

rve = Rve(dim_x=1, dim_y=1, dim_z=1)

elem = Ellipsoid(a_x=0.15, a_y=0.31, a_z=0.4)
elli = elem.generate()

raster = rasterPhase(phase=Phase(shape=elli), rve=rve, grid=[5, 5, 5])

compound = cq.Compound.makeCompound(
    [solid for phase in raster for solid in phase.solids]
)
cq.exporters.export(compound, "compound.step")

mesh(
    mesh_file="compound.step",
    listPhases=raster,
    size=0.03,
    order=1,
    output_file="rasterEllipsoid.vtk",
)
