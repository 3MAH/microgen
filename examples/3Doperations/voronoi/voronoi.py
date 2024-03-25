from pathlib import Path

import cadquery as cq

from microgen import Neper, Phase, mesh

# # We import the Polyhedra from Neper tessellation file
# listPolyhedra, seed, vertices, edges, faces, polys = parseNeper("test1")

# Revel = Rve(dim=1)
# phases = []

# for polyhedron in listPolyhedra:
#     elem = Polyhedron(
#         center=(
#             polyhedron["original"][0],
#             polyhedron["original"][1],
#             polyhedron["original"][2],
#         ),
#         dic=polyhedron,
#     )
#     phases.append(Phase(shape=elem.generate()))

# compound = cq.Compound.makeCompound([phase.shape for phase in phases])
# cq.exporters.export(compound, "compound.step")

# mesh(
#     mesh_file="compound.step",
#     listPhases=phases,
#     size=0.05,
#     order=1,
#     output_file="Voronoi.vtk",
# )

tess_file = str(Path(__file__).parent / "test1.tess")
polyhedra = Neper.generateVoronoiFromTessFile(tess_file)

shapes = [poly.generate() for poly in polyhedra]

compound = cq.Compound.makeCompound(shapes)
step_file = str(Path(__file__).parent / "compound.step")
cq.exporters.export(compound, step_file)

phases = [Phase(shape=shape) for shape in shapes]

vtk_file = str(Path(__file__).parent / "Voronoi.vtk")
mesh(
    mesh_file=step_file,
    listPhases=phases,
    size=0.05,
    order=1,
    output_file=vtk_file,
)
