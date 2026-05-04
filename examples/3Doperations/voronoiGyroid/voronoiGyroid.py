from pathlib import Path

from microgen import Neper, Phase, Tpms, mesh
from microgen.cad import make_compound
from microgen.mesh import MeshOptions
from microgen.shape import surface_functions

# We import the Polyhedra from Neper tessellation file
tess_file = str(Path(__file__).parent / "test1.tess")
polyhedra = Neper.generateVoronoiFromTessFile(tess_file)

gyroid = Tpms(
    center=(0.5, 0.5, 0.5),
    surface_function=surface_functions.gyroid,
    offset=0.2,
)
gyroid = gyroid.generate(type_part="sheet").translate((0.5, 0.5, 0.5))

phases = []
for polyhedron in polyhedra:
    shape = polyhedron.generate()
    phases.append(Phase(shape=shape.intersect(gyroid)))

compound = make_compound([phase.shape for phase in phases])
step_file = str(Path(__file__).parent / "compound.step")
compound.export_step(step_file)

vtk_file = str(Path(__file__).parent / "Gyroid-voro.vtk")
mesh(
    mesh_file=step_file,
    list_phases=phases,
    options=MeshOptions(size=0.05, order=1, output_file=vtk_file),
)
