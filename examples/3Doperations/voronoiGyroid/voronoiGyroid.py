from pathlib import Path

from microgen import Neper, Phase, Tpms, mesh
from microgen.cad import make_compound
from microgen.shape import surface_functions

# We import the Polyhedra from Neper tessellation file
tess_file = str(Path(__file__).parent / "test1.tess")
polyhedra = Neper.voronoi_from_tess_file(tess_file)

gyroid = Tpms(
    center=(0.5, 0.5, 0.5),
    surface_function=surface_functions.gyroid,
    offset=0.2,
)
gyroid = gyroid.generate_cad(type_part="sheet").translate((0.5, 0.5, 0.5))

phases = []
for polyhedron in polyhedra:
    shape = polyhedron.generate_cad()
    phases.append(Phase.from_cad(shape.intersect(gyroid)))

compound = make_compound([phase.cad for phase in phases])
step_file = str(Path(__file__).parent / "compound.step")
compound.export_step(step_file)

vtk_file = str(Path(__file__).parent / "Gyroid-voro.vtk")
mesh(
    mesh_file=step_file,
    list_phases=phases,
    size=0.05,
    order=1,
    output_file=vtk_file,
)
