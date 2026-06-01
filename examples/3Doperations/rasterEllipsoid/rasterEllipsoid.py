from pathlib import Path

from microgen import Ellipsoid, Phase, Rve, mesh, raster_phase
from microgen.cad import make_compound_from_solids

rve = Rve(dim=1)

elem = Ellipsoid(radii=(0.15, 0.31, 0.4))
elli = elem.generate_cad()

raster = raster_phase(phase=Phase.from_cad(elli), rve=rve, grid=[5, 5, 5])

compound = make_compound_from_solids(
    [solid.wrapped for phase in raster for solid in phase.cad.solids()]
)
step_file = str(Path(__file__).parent / "compound.step")
compound.export_step(step_file)

vtk_file = str(Path(__file__).parent / "rasterEllipsoid.vtk")
mesh(
    mesh_file=step_file,
    list_phases=raster,
    size=0.03,
    order=1,
    output_file=vtk_file,
)
