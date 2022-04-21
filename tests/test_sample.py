import microgen
import cadquery as cq
import meshio
import numpy as np

import pytest

def test_tpms():
    rve = microgen.Rve(dim_x=1, dim_y=1, dim_z=1, size_mesh=0.03)
    bounding_sphere_radius = 1.1*np.sqrt((0.5*rve.dx)**2 + (0.5*rve.dy)**2 + (0.5*rve.dz)**2)

    function = "cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)"
    height = 0.3

    tpms = microgen.Custom(rve, function, height)
    assert tpms.eval([0, 0, 0]) == 3 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    
    max_edges_size_at_feature_edges = 4
    tpms = microgen.Hyperboloid(max_edges_size_at_feature_edges)
    tpms.getFeatures()
    assert tpms.eval([0, 0, 0]) == -0.5**2
    assert tpms.get_bounding_sphere_squared_radius() == (1 + 0.5)**2 + 1**2
    tpms = microgen.SchwarzP(rve, height)
    assert tpms.eval([0, 0, 0]) == 3 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.SchwarzD(rve, height)
    assert tpms.eval([0, 0, 0]) == 0 + 0 + 0 + 0 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.Neovius(rve, height)
    assert tpms.eval([0, 0, 0]) == (3 + 1 + 1) + (4*1*1*1) + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.SchoenIWP(rve, height)
    assert tpms.eval([0, 0, 0]) == 2*(1+1+1) - (1+1+1) + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.SchoenFRD(rve, height)
    assert tpms.eval([0, 0, 0]) == 4 - (1+1+1) + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.FischerKochS(rve, height)
    assert tpms.eval([0, 0, 0]) == 0 + 0 + 0 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.PMY(rve, height)
    assert tpms.eval([0, 0, 0]) == 2 + 0 + 0 + 0 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.HoneyComb(rve, height)
    assert tpms.eval([0, 0, 0]) == 0 + 0 + 1 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    tpms = microgen.Gyroid(rve, height)
    assert tpms.eval([0, 0, 0]) == 1
    # assert tpms.eval([1, 0, 0]) == 0 + 0 + 0 + height
    assert tpms.get_bounding_sphere_squared_radius() == bounding_sphere_radius


def test_shapes():
    rve = microgen.Rve(a=1, b=1, c=1, size_mesh=0.03)

    elem = microgen.BasicGeometry(number=0, shape='ellipsoid',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"a_x": 0.15,
                                              "a_y": 0.31,
                                              "a_z": 0.4})
    shape = elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='sphere',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"radius": 0.15})
    shape = elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='box',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"dim_x": 0.15,
                                              "dim_y": 0.31,
                                              "dim_z": 0.4})
    shape = elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='capsule',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"height": 0.5,
                                              "radius": 0.1})
    shape = elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='cylinder',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"height": 0.5,
                                              "radius": 0.1})
    shape = elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='extrudedpolygon',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"listCorners": [(0, 0), (0, 1), (1, 1), (1, 0)],
                                              "height": 0.3})
    shape = elem.generate()

   dic = {'original': [0.107874084791, 0.618217780057, 0.938426948697], 
           'faces': [{'vertices': [0, 1, 2, 3, 4]}, 
                     {'vertices': [5, 6, 2, 1, 7, 8]}, 
                     {'vertices': [9, 5, 6, 10, 11]}, 
                     {'vertices': [7, 1, 0, 12]}, 
                     {'vertices': [9, 11, 13, 14, 15]}, 
                     {'vertices': [13, 4, 0, 12, 14]}, 
                     {'vertices': [5, 8, 15, 9]}, 
                     {'vertices': [2, 3, 10, 6]}, 
                     {'vertices': [4, 3, 10, 11, 13]}, 
                     {'vertices': [15, 8, 7, 12, 14]}], 
                      'vertices': [[-0.0, 0.523214488307, 0.820162121852], 
                                   [-0.0, 0.513287784699, 1.0], 
                                   [-0.0, 0.740463797741, 1.0], 
                                   [-0.0, 0.755067696893, 0.869556881388], 
                                   [-0.0, 0.648014540172, 0.772167416819], 
                                   [0.251073490918, 0.663754005059, 1.0], 
                                   [0.184284417298, 0.738145126407, 1.0], 
                                   [0.151624118287, 0.48362847958, 1.0], 
                                   [0.255188032397, 0.570803855492, 1.0], 
                                   [0.251051495806, 0.635892401699, 0.86011055724], 
                                   [0.143585076339, 0.754776197774, 0.856023965986], 
                                   [0.185032400205, 0.697452142104, 0.800044686732], 
                                   [0.152116967029, 0.493611691791, 0.817391823705], 
                                   [0.142129671376, 0.647943344513, 0.758969402247], 
                                   [0.190371122851, 0.526664864561, 0.80112086188], 
                                   [0.253676937949, 0.57680436693, 0.861207628252]]}
    elem = microgen.Phase.BasicGeometry(number=0, shape='polyhedron',
                                        xc=0, yc=0, zc=0,
                                        psi=0, theta=0, phi=0,
                                        param_geom={"dic": dic})
    shape = elem.generate()



    elem = microgen.Phase.BasicGeometry(number=0, shape='fake',
                                  xc=0.5, yc=0.5, zc=0.5,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"geom": 0})
    pytest.raises(ValueError, elem.generate)

    elem = microgen.BasicGeometry(number=0, shape='tpms',
                                  xc=0.5, yc=0.5, zc=0.5,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"type_surface": "custom",
                                              "type_part": "skeletal",
                                              "thickness": 0.1,
                                              "function": 'cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)'},
                                  path_data='tests/data')
    elem.geometry.createSurfaces(rve=rve,
                                 sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                                 path_data='tests/data')
    shape = elem.generate(rve=rve)

    elem = microgen.BasicGeometry(number=0, shape='tpms',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"type_surface": "gyroid",
                                              "type_part": "sheet",
                                              "thickness": 0.3},
                                  path_data='tests/data')
    elem.geometry.createSurfaces(rve=rve,
                                 sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                                 path_data='tests/data')
    shape = elem.generate(rve=rve)

    raster = microgen.rasterShapeList(cqShapeList=[shape], rve=rve, grid=[5, 5, 5])

    compound = cq.Compound.makeCompound(raster[0])
    cq.exporters.export(compound, 'tests/data/compound.step')

    microgen.mesh(mesh_file='tests/data/compound.step', listPhases=raster[1], size=0.03, order=1)
    

def test_octettruss():

    # fichier
    NPhases_file = "examples/octetTruss/test_octet.dat"
    microgen.removeEmptyLines(NPhases_file)

    dt = np.dtype([('number', int), ('shape', np.str_, 10),
                ('xc', np.float64), ('yc', np.float64), ('zc', np.float64),
                ('psi', np.float64), ('theta', np.float64), ('phi', np.float64),
                ('height', np.float64), ('radius', np.float64)])
    # précision du type des données
    number, shape, xc, yc, zc, psi, theta, phi, height, radius, = np.loadtxt(NPhases_file, dtype=dt,
                                                                             usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                                                                             skiprows=1, unpack=True, ndmin=1)

    # Size of the mesh
    size_mesh = 0.03
    a = 1.0
    b = 1.0
    c = 1.0

    revel = microgen.Rve(a, b, c, size_mesh)
    listPhases = []
    listPeriodicPhases = []
    n = len(xc)

    for i in range(0, n):
        elem = microgen.BasicGeometry(number=number[i], shape=shape[i],
                            xc=xc[i], yc=yc[i], zc=zc[i],
                            psi=psi[i], theta=theta[i], phi=phi[i],
                            param_geom={"height": height[i],
                                        "radius": radius[i]},
                            path_data='')
        listPhases.append(elem.generate())

    for phase_elem in listPhases:
        print(phase_elem)
        periodicPhase = microgen.periodic(cqshape=phase_elem, rve=revel)
        listPeriodicPhases.append(periodicPhase)

    phases_cut = microgen.cutParts(cqShapeList=[s[0] for s in listPeriodicPhases], reverseOrder=False)
    compound = cq.Compound.makeCompound(phases_cut[0])

    cq.exporters.export(compound, 'tests/data/compound.step')
    microgen.meshPeriodic(mesh_file='tests/data/compound.step', rve=revel, listPhases=phases_cut[1], size=0.03, order=1)


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
