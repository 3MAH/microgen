import cadquery as cq

from microgen import Neper, Phase, Tpms, mesh
from microgen.shape import surface_functions

# We import the Polyhedra from Neper tessellation file
polyhedra = Neper.generateVoronoiFromTessFile("test1")

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
cq.exporters.export(compound, "compound.step")

mesh(
    mesh_file="compound.step",
    listPhases=phases,
    size=0.05,
    order=1,
    output_file="Gyroid-voro.vtk",
)
