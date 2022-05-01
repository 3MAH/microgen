import microgen
import cadquery as cq
import numpy as np


def test_misc():
    assert microgen.is_function_valid("cos(2*pi*x) + 3*sin(4*pi*x)") is True
    assert microgen.is_function_valid("co(2*pi*x) + 3*sin(4*pi*x)") is False
    assert microgen.is_function_valid("cos(2pi*x) + 3*sin(4*pi*x)") is False
    
    listPolyhedra, seed, vertices, edges, faces, polys = microgen.parseNeper('examples/Voronoi/test1')

def test_cut():

    
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



    # FUNCTIONS
    # lanceNeper

    # MATERIAL
    # readSections
