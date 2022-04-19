import microgen

import numpy as np

def test_tpms():
    rve = microgen.Rve(a=1, b=1, c=1, size_mesh=0.03)
    function = "cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)"
    height = 0.3

    tpms = microgen.Custom(rve, function, height)
    assert tpms.eval([0, 0, 0]) == 3 + height
    assert tpms.get_bounding_sphere_squared_radius() == 1.1 * np.sqrt((0.5**2 + 0.5**2 + 0.5**2))
    
    max_edges_size_at_feature_edges = 4
    tpms = microgen.Hyperboloid(max_edges_size_at_feature_edges)
    assert tpms.eval([0, 0, 0]) == -0.5**2
    tpms = microgen.SchwarzP(rve, height)
    assert tpms.eval([0, 0, 0]) == 3 + height
    tpms = microgen.SchwarzD(rve, height)
    assert tpms.eval([0, 0, 0]) == 0 + 0 + 0 + 0 + height
    tpms = microgen.Neovius(rve, height)
    assert tpms.eval([0, 0, 0]) == (3 + 1 + 1) + (4*1*1*1) + height
    tpms = microgen.SchoenIWP(rve, height)
    assert tpms.eval([0, 0, 0]) == 2*(1+1+1) - (1+1+1) + height
    tpms = microgen.SchoenFRD(rve, height)
    assert tpms.eval([0, 0, 0]) == 4 - (1+1+1) + height
    tpms = microgen.FischerKochS(rve, height)
    assert tpms.eval([0, 0, 0]) == 0 + 0 + 0 + height
    tpms = microgen.PMY(rve, height)
    assert tpms.eval([0, 0, 0]) == 2 + 0 + 0 + 0 + height
    tpms = microgen.HoneyComb(rve, height)
    assert tpms.eval([0, 0, 0]) == 0 + 0 + 1 + height
    tpms = microgen.Gyroid(rve, height)
    assert tpms.eval([0, 0, 0]) == 1
    # assert tpms.eval([1, 0, 0]) == 0 + 0 + 0 + height


    assert microgen.is_function_valid("cos(2*pi*x) + 3*sin(4*pi*x)") == True
    assert microgen.is_function_valid("co(2*pi*x) + 3*sin(4*pi*x)") == False
    assert microgen.is_function_valid("cos(2pi*x) + 3*sin(4*pi*x)") == False

    microgen.BasicGeometry(number=0, shape='ellipsoid',
                           xc=0, yc=0, zc=0,
                           psi=0, theta=0, phi=0,
                           param_geom={"a1": 0.15,
                                       "a2": 0.31,
                                       "a3": 0.4},
                           path_data='')

    microgen.BasicGeometry(number=0, shape='sphere',
                           xc=0, yc=0, zc=0,
                           psi=0, theta=0, phi=0,
                           param_geom={"radius": 0.15},
                           path_data='')

    microgen.BasicGeometry(number=0, shape='box',
                           xc=0, yc=0, zc=0,
                           psi=0, theta=0, phi=0,
                           param_geom={"a1": 0.15,
                                       "a2": 0.31,
                                       "a3": 0.4},
                           path_data='')

    microgen.BasicGeometry(number=0, shape='capsule',
                           xc=0, yc=0, zc=0,
                           psi=0, theta=0, phi=0,
                           param_geom={"height": 0.5,
                                       "radius": 0.1},
                           path_data='')

    microgen.BasicGeometry(number=0, shape='cylinder',
                           xc=0, yc=0, zc=0,
                           psi=0, theta=0, phi=0,
                           param_geom={"height": 0.5,
                                       "radius": 0.1},
                           path_data='')

    microgen.BasicGeometry(number=0, shape='extrudedpolygon',
                           xc=0, yc=0, zc=0,
                           psi=0, theta=0, phi=0,
                           param_geom={"listCorners": [(0, 0), (0, 1), (1, 1), (1, 0)],
                                       "height": 0.3},
                           path_data='')

    # microgen.BasicGeometry(number=0, shape='polyhedron',
    #                        xc=0, yc=0, zc=0,
    #                        psi=0, theta=0, phi=0,
    #                        param_geom={},
    #                        path_data='')

    microgen.BasicGeometry(number=0, shape='tpms',
                           xc=0, yc=0, zc=0,
                           psi=0, theta=0, phi=0,
                           param_geom={"type_surface": "custom",
                                       "type_part": "sheet",
                                       "thickness": 0.3,
                                       "function": 'cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)'},
                           path_data='')

    microgen.BasicGeometry(number=0, shape='tpms',
                           xc=0, yc=0, zc=0,
                           psi=0, theta=0, phi=0,
                           param_geom={"type_surface": "gyroid",
                                       "type_part": "sheet",
                                       "thickness": 0.3},
                           path_data='')