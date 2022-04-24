import cadquery as cq
from microgen import Rve, fuseParts, BasicGeometry

def repeatGeometry(unit_geom, rve, grid):
    import time

    start = time.time()

    xyz_repeat = cq.Assembly()
    for i_x in range(grid["x"]):
        for i_y in range(grid["y"]):
            for i_z in range(grid["z"]):
                xyz_repeat.add(unit_geom, 
                               loc=cq.Location(cq.Vector(i_x*rve.dim_x, 
                                                         i_y*rve.dim_y, 
                                                         i_z*rve.dim_z)))


    end = time.time()

    print(end - start)
    
    return xyz_repeat.toCompound()

# Size of the mesh
size_mesh = 0.03
a = 1.0
b = 1.0
c = 1.0

rve = Rve(dim_x=a, dim_y=b, dim_z=c, size_mesh=size_mesh)

# unit_geom = cq.importers.importStep('../../octetTruss/compound.step')
unit_geom = cq.importers.importStep('../../triplyPeriodicMinimalSurfaces/sheet.step')
unit_geom = cq.Shape(unit_geom.val().wrapped)

final_geom = repeatGeometry(unit_geom, rve, grid={"x": 3, 
                                                  "y": 3, 
                                                  "z": 3})

cq.exporters.export(final_geom, 'compound.stl')