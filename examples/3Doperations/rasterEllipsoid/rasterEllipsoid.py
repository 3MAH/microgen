from pathlib import Path

from microgen import Ellipsoid, Phase, Rve, mesh, raster_phase
from microgen.cad import make_compound_from_solids
from microgen.mesh import MeshOptions

rve = Rve(dim=1)

elem = Ellipsoid(radii=(0.15, 0.31, 0.4))
elli = elem.generate()

raster = raster_phase(phase=Phase(shape=elli), rve=rve, grid=[5, 5, 5])

compound = make_compound_from_solids(
    [solid for phase in raster for solid in phase.solids]
)
step_file = str(Path(__file__).parent / "compound.step")
compound.export_step(step_file)

vtk_file = str(Path(__file__).parent / "rasterEllipsoid.vtk")
mesh(
    mesh_file=step_file,
    list_phases=raster,
    options=MeshOptions(size=0.03, order=1, output_file=vtk_file),
)
