import cadquery as cq
from microgen import parseNeper, Rve, BasicGeometry, mesh

# We import the Polyhedra from Nepoer tesselation file
listPolyhedra, seed, vertices, edges, faces, polys = parseNeper('test1')

# We fix the size of the mesh and dimensions of the RVE
size_mesh = 0.05
a = 1.0
b = 1.0
c = 1.0

periodicity = 0

Revel = Rve(a, b, c)
phases = []
listSolids = []

for polyhedron in listPolyhedra:
    elem = BasicGeometry(number=0, shape='polyhedron',
                         xc=polyhedron["original"][0],
                         yc=polyhedron["original"][1],
                         zc=polyhedron["original"][2],
                         psi=0., theta=0., phi=0.,
                         param_geom={"dic": polyhedron}, path_data='.')
    temp = elem.generate()
    phases.append(temp)

compound = cq.Compound.makeCompound(phases)
cq.exporters.export(compound, 'compound.step')

occ_solids_list = [s.Solids() for s in phases]
mesh(mesh_file='compound.step', listPhases=occ_solids_list, size=size_mesh, order=1, output_file='Mesh.vtk')
