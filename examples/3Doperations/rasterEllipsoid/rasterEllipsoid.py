from pathlib import Path

import cadquery as cq

from microgen import Ellipsoid, Phase, Rve, mesh, rasterPhase

rve = Rve(dim=1)

elem = Ellipsoid(a_x=0.15, a_y=0.31, a_z=0.4)
elli = elem.generate()

raster = rasterPhase(phase=Phase(shape=elli), rve=rve, grid=[5, 5, 5])

compound = cq.Compound.makeCompound(
    [solid for phase in raster for solid in phase.solids]
)
step_file = str(Path(__file__).parent / "compound.step")
cq.exporters.export(compound, step_file)

vtk_file = str(Path(__file__).parent / "rasterEllipsoid.vtk")
mesh(
    mesh_file=step_file,
    listPhases=raster,
    size=0.03,
    order=1,
    output_file=vtk_file,
)
