import microgen
import cadquery as cq
import meshio

import numpy as np

def test_tpms():
    rve = microgen.Rve(a=1, b=1, c=1, size_mesh=0.03)
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

    # elem = microgen.BasicGeometry(number=0, shape='polyhedron',
    #                               xc=0, yc=0, zc=0,
    #                               psi=0, theta=0, phi=0,
    #                               param_geom={})
    # shape = elem.generate()

    # elem = microgen.BasicGeometry(number=0, shape='tpms',
    #                               xc=0.5, yc=0.5, zc=0.5,
    #                               psi=0, theta=0, phi=0,
    #                               param_geom={"type_surface": "custom",
    #                                           "type_part": "skeletal",
    #                                           "thickness": 0.1,
    #                                           "function": 'cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)'},
    #                               path_data='tests/data')
    # elem.geometry.createSurfaces(rve=rve,
    #                              sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
    #                              path_data='tests/data')
    # shape = elem.generate(rve=rve)

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
    assert microgen.is_function_valid("cos(2*pi*x) + 3*sin(4*pi*x)") == True
    assert microgen.is_function_valid("co(2*pi*x) + 3*sin(4*pi*x)") == False
    assert microgen.is_function_valid("cos(2pi*x) + 3*sin(4*pi*x)") == False
    
    listPolyhedra, seed, vertices, edges, faces, polys = microgen.parseNeper('examples/Voronoi/test1')


    # FUNCTIONS
    # lanceNeper
    # cutPhaseByShapeList
    # cutPhasesByShape

    # MATERIAL
    # readSections