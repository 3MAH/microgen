import cadquery as cq
from microgen import parseNeper, Tpms, Polyhedron, Phase, mesh
import microgen

# We import the Polyhedra from Nepoer tesselation file
listPolyhedra = parseNeper("test1")[0]

gyroid = Tpms(
    center=(0.5, 0.5, 0.5),
    surface_function=microgen.shape.tpms.gyroid,
    type_part="sheet",
    thickness=0.2,
    path_data="data",
)
gyroid = gyroid.generate(sizeMesh=0.03, minFacetAngle=20.0, maxRadius=0.03).translate(
    (0.5, 0.5, 0.5)
)

phases = []
i = 0
for polyhedron in listPolyhedra:
    elem = Polyhedron(
        center=(
            polyhedron["original"][0],
            polyhedron["original"][1],
            polyhedron["original"][2],
        ),
        dic=polyhedron,
    )
    elem = elem.generate()
    phases.append(Phase(shape=elem.intersect(gyroid)))

compound = cq.Compound.makeCompound([phase.shape for phase in phases])
cq.exporters.export(compound, "compound.step")

mesh(
    mesh_file="compound.step",
    listPhases=phases,
    size=0.05,
    order=1,
    output_file="Gyroid-voro.vtk",
)
