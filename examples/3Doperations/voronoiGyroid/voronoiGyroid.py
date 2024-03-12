from pathlib import Path

import cadquery as cq

from microgen import Neper, Phase, Tpms, mesh
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

compound = cq.Compound.makeCompound([phase.shape for phase in phases])
step_file = str(Path(__file__).parent / "compound.step")
cq.exporters.export(compound, step_file)

vtk_file = str(Path(__file__).parent / "Gyroid-voro.vtk")
mesh(
    mesh_file=step_file,
    listPhases=phases,
    size=0.05,
    order=1,
    output_file=vtk_file,
)
