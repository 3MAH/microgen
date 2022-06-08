import cadquery as cq
from microgen import parseNeper, Rve, Polyhedron, mesh, Phase

# We import the Polyhedra from Nepoer tesselation file
listPolyhedra, seed, vertices, edges, faces, polys = parseNeper('test1')

Revel = Rve(dim_x=1, dim_y=1, dim_z=1)
phases = []
listSolids = []

for polyhedron in listPolyhedra:
    elem = Polyhedron(center=(polyhedron["original"][0], polyhedron["original"][1], polyhedron["original"][2]),
                      dic=polyhedron)
    phases.append(Phase(shape=elem.generate()))

compound = cq.Compound.makeCompound([phase.shape for phase in phases])
cq.exporters.export(compound, 'compound.step')

mesh(mesh_file='compound.step', listPhases=phases, size=0.05, order=1, output_file='Voronoi.vtk')
