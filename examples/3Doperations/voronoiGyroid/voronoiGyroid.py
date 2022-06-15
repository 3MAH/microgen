import cadquery as cq
from microgen import parseNeper, Tpms, Polyhedron, Phase, mesh, Neper
from microgen.shape import tpms

# We import the Polyhedra from Nepoer tesselation file
polyhedra = Neper.generateVoronoiFromTessFile("test1")

gyroid = Tpms(
    center=(0.5, 0.5, 0.5),
    surface_function=tpms.gyroid,
    type_part="sheet",
    thickness=0.2,
    path_data="data",
)
gyroid = gyroid.generate().translate(
    (0.5, 0.5, 0.5)
)

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
