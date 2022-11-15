from microgen import Tpms, tpms, fuseShapes

import numpy as np
import pyvista as pv
import cadquery as cq

import sys

import itertools

def createShell(mesh: pv.PolyData) -> cq.Shell:
    list_of_triangles = mesh.faces.reshape(-1, 4)[:, 1:]
    list_of_triangles = np.c_[list_of_triangles, list_of_triangles[:, 0]]

    faces = []
    for i in range(len(list_of_triangles)):
        ixs = list_of_triangles[i]
        print(f"\t{i}/{len(list_of_triangles)}", end="\r")
        lines = []
        for v1, v2 in zip(ixs[:], ixs[1:]):
            vertice_coords1 = mesh.points[v1]
            vertice_coords2 = mesh.points[v2]
            lines.append(
                cq.Edge.makeLine(
                    cq.Vector(*vertice_coords1), cq.Vector(*vertice_coords2)
                )
            )
        wire = cq.Wire.assembleEdges(lines)
        faces.append(cq.Face.makeFromWires(wire))
    return cq.Shell.makeShell(faces)

def h(x, y, z, cell_size):
    hmin = 0.05
    hmax = 0.25
    L = cell_size[0]

    return (hmax - hmin) * x / L + 0.5 * (hmin + hmax)

cell_size = (5, 5, 1)
repeat_cell = (5, 5, 1)
resolution = 20

res = (resolution * repeat_cell[0], resolution * repeat_cell[1], resolution * repeat_cell[2])
x_min, y_min, z_min = -0.5 * cell_size[0], -0.5 * cell_size[1], -0.5 * cell_size[2]
grid = pv.UniformGrid(
    dims=res,
    spacing=(cell_size[0] / (res[0] - 1),
             cell_size[1] / (res[1] - 1),
             cell_size[2] / (res[2] - 1)),
    origin=(x_min, y_min, z_min),
)
x, y, z = grid.points.T

surface = tpms.gyroid(x, y, z, repeat_cell, cell_size)
# surface = tpms.schwarzP(x, y, z, repeat_cell, cell_size)
# surface = tpms.schwarzD(x, y, z, repeat_cell, cell_size)
# surface = tpms.schoenFRD(x, y, z, repeat_cell, cell_size)
# surface = tpms.schoenIWP(x, y, z, repeat_cell, cell_size)
# surface = tpms.neovius(x, y, z, repeat_cell, cell_size) # Error
# surface = tpms.honeycomb(x, y, z, repeat_cell, cell_size) # bizarre
# surface = tpms.pmy(x, y, z, repeat_cell, cell_size) # bizarre
# surface = tpms.fischerKochS(x, y, z, repeat_cell, cell_size) # bizarre
# surface = tpms.lidinoid(x, y, z, repeat_cell, cell_size) # ERROR

thickness = 2 * np.pi * h(x, y, z, cell_size)

surf_p = surface + thickness / 2.
mesh_p = grid.contour(isosurfaces=[0], scalars=surf_p)
mesh_p.smooth(n_iter=100)
mesh_p.clean(inplace=True)
shell_p = createShell(mesh_p)
print("shell_p")

surf_m = surface - thickness / 2.
mesh_m = grid.contour(isosurfaces=[0], scalars=surf_m)
mesh_m.smooth(n_iter=100)
mesh_m.clean(inplace=True)
shell_m = createShell(mesh_m)
print("shell_m")


box = cq.Workplane("front").box(cell_size[0], cell_size[1], cell_size[2])
box_cut = box.split(shell_p).split(shell_m)

print("cut")

solids = box_cut.solids().all()
print(len(solids))

if len(solids) % 2 == 0: # assume that the number of solids is odd
    print("ERROR")
    sys.exit()

middle_solid_number = (len(solids) - 1) // 2
if (middle_solid_number) % 2 == 0: 
    sheet = [cq.Shape(solids[even].val().wrapped) for even in range(0, len(solids), 2)]
    skeletal = [cq.Shape(solids[odd].val().wrapped) for odd in range(1, len(solids), 2)]
else:
    skeletal = [cq.Shape(solids[even].val().wrapped) for even in range(0, len(solids), 2)]
    sheet = [cq.Shape(solids[odd].val().wrapped) for odd in range(1, len(solids), 2)]

print("sheet and skeletal separated")

sheet = cq.Compound.makeCompound(sheet)
sheet.exportStl("sheet.stl")
skeletal = [cq.Compound.makeCompound(skeletal[:len(skeletal)//2]), cq.Compound.makeCompound(skeletal[len(skeletal)//2:])]

print("converting to polydata")
pd_sheet = pv.PolyData(sheet.toVtkPolyData(tolerance=0.01, angularTolerance=0.1, normals=True))
pd_skeletal_m = pv.PolyData(skeletal[0].toVtkPolyData(tolerance=0.01, angularTolerance=0.1, normals=True))
pd_skeletal_p = pv.PolyData(skeletal[1].toVtkPolyData(tolerance=0.01, angularTolerance=0.1, normals=True))

p = pv.Plotter()
p.add_mesh(pd_sheet, color='w')
p.add_mesh(pd_skeletal_m, color='b')
p.add_mesh(pd_skeletal_p, color='r')
p.camera_position = 'xy'
p.parallel_projection = True 
p.show_axes()
p.show()