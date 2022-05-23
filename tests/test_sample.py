import microgen
import cadquery as cq
import numpy as np


def test_misc():    
    listPolyhedra, seed, vertices, edges, faces, polys = microgen.parseNeper('examples/Voronoi/test1')

    microgen.removeEmptyLines("fake_file.txt")

    microgen.MatSection(number=0,
                        name="test",
                        umat_name="test",
                        psi_mat=0,
                        theta_mat=0,
                        phi_mat=0,
                        nprops=0,
                        nstatev=0,
                        props=[])

def test_operations():

    
    elem = microgen.BasicGeometry(number=0, shape='box',
                                  xc=0.5, yc=0.5, zc=0.5,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"dim_x": 1,
                                              "dim_y": 1,
                                              "dim_z": 1})
    shape1 = elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='box',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"dim_x": 0.5,
                                              "dim_y": 0.5,
                                              "dim_z": 0.5})
    shape2 = elem.generate()

    microgen.cutPhaseByShapeList(shape1, [shape2])

    microgen.cutPhasesByShape([shape1], shape2)

    rve = microgen.Rve(dim_x=1, dim_y=1, dim_z=1, size_mesh=0.03)
    microgen.repeatGeometry(shape1, rve, grid={"x": 2, "y": 2, "z": 2})



    # EXTERNAL
    # lanceNeper

    # MATERIAL
    # readSections
