from microgen import tpms, repeatShape, Rve

import numpy as np
import pyvista as pv

def offset(x, y, z):
    hmin = 0.05
    hmax = 0.45
    length = x[-1] - x[0]
    return (hmax - hmin) * x / length + 0.5 * (hmin + hmax)

surface = tpms.gyroid
# surface = tpms.schwarzP
# surface = tpms.schwarzD
# surface = tpms.schoenFRD
# surface = tpms.schoenIWP
# surface = tpms.neovius
# surface = tpms.honeycomb # bizarre
# surface = tpms.pmy # bizarre
# surface = tpms.fischerKochS # bizarre
# surface = tpms.lidinoid # ERROR


cell_size = np.array([5, 1, 1])
repeat_cell = np.array([5, 1, 1])
resolution = 20

preview = False

geometry = tpms.Tpms(surface_function=surface,
                     offset=offset, 
                     repeat_cell=repeat_cell, 
                     cell_size=cell_size)
sheet = geometry.generate(type_part=tpms.SHEET, 
                          resolution=resolution, 
                          verbose=True, 
                          preview=preview)


sheet = repeatShape(unit_geom=sheet, rve=Rve(5, 1, 1), grid=(1, 5, 5))
sheet.exportStl("sheet.stl")

# print("converting to polydata")
pd_sheet = pv.PolyData(sheet.toVtkPolyData(tolerance=0.01, angularTolerance=0.1, normals=True))
# pd_skeletals = pv.PolyData(skeletals.toVtkPolyData(tolerance=0.01, angularTolerance=0.1, normals=True))

p = pv.Plotter()
p.add_mesh(pd_sheet, color='w')
# p.add_mesh(pd_skeletals, color='b')
p.camera_position = 'xy'
p.parallel_projection = True 
p.show_axes()
p.show()