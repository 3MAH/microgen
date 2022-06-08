import cadquery as cq
from microgen import Tpms, tpms
import numpy as np

surfaces = [tpms.neovius, tpms.gyroid, tpms.schwarzD,
            tpms.schwarzP, tpms.honeycomb, tpms.schoenIWP,
            tpms.schoenFRD, tpms.fischerKochS, tpms.pmy]

assembly = cq.Assembly()

i = 0
n_col = 3
n_row = np.ceil(len(surfaces) / n_col)
for surface in surfaces:
    i_x = i % n_col
    i_y = i // n_col

    elem = Tpms(center=(0.5, 0.5, 0.5),
                orientation=(90, 90, 90),
                surface_function=surface,
                type_part='sheet',
                thickness=0.,
                path_data="data")
    elem.createSurface(sizeMesh=0.03, minFacetAngle=10., maxRadius=0.03)
    center = (1.2 * (i_x - 0.5 * (n_col - 1)), -1.2 * (i_y - 0.5 * (n_row - 1)), 0)
    assembly.add(elem.generateSurface(), loc=cq.Location(cq.Vector(center)))

    i = i + 1
    print(" {}/{}".format(i, len(surfaces)), end='\r')

cq.exporters.export(assembly.toCompound(), 'tpms.stl')
