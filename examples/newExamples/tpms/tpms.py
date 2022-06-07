import cadquery as cq
from microgen import BasicGeometry, Rve
import time
import numpy as np

total_start = time.time()

surfaces = ["neovius", "gyroid", "schwarzD", 
            "schwarzP", "honeycomb", "schoenIWP",
            "schoenFRD", "fischerKochS", "pmy"]

rve = Rve(dim_x=1, dim_y=1, dim_z=1)

assembly = cq.Assembly()

i = 0
n_col = 3
n_row = np.ceil(len(surfaces)/n_col)
for surface in surfaces:
    i_x = i%n_col
    i_y = i//n_col

    start = time.time()

    elem = BasicGeometry(number=i, shape='tpms',
                         xc=0.5, yc=0.5, zc=0.5,
                         psi=90, theta=90, phi=90,
                         param_geom={"type_surface": surface,
                                     "type_part": 'sheet',
                                     "thickness": 0.1},
                         path_data="data")
    elem.geometry.createSurfaces(rve=rve,
                                 sizeMesh=0.05, minFacetAngle=20., maxRadius=0.05,
                                 path_data='data', verbose=False)
    assembly.add(elem.generate(rve=rve), 
                 loc=cq.Location(cq.Vector( 1.2*(i_x - 0.5*(n_col - 1)), 
                                           -1.2*(i_y - 0.5*(n_row - 1)), 
                                            0
                                          )
                                )
                )

    end = time.time()

    print(surface, end - start)

    i = i + 1

total_end = time.time()

print("Total time : {} s".format(total_end - total_start))

# show_object(compound) # run with cq-editor

cq.exporters.export(assembly.toCompound(), 'tpms.stl')
# cq.exporters.export(assembly.toCompound(), 'tpms.step')