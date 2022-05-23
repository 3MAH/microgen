import microgen
import cadquery as cq
import pytest
import numpy as np

def test_shapes():
    rve = microgen.Rve(dim_x=1, dim_y=1, dim_z=1, size_mesh=0.03)

    elem = microgen.BasicGeometry(number=0, shape='ellipsoid',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"a_x": 0.15,
                                              "a_y": 0.31,
                                              "a_z": 0.4})
    ellipsoid = elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='sphere',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"radius": 0.15})
    elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='box',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"dim_x": 0.15,
                                              "dim_y": 0.31,
                                              "dim_z": 0.4})
    elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='capsule',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"height": 0.5,
                                              "radius": 0.1})
    elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='cylinder',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"height": 0.5,
                                              "radius": 0.1})
    elem.generate()

    elem = microgen.BasicGeometry(number=0, shape='extrudedpolygon',
                                  xc=0, yc=0, zc=0,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"listCorners": [(0, 0), (0, 1), (1, 1), (1, 0)],
                                              "height": 0.3})
    elem.generate()

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
    elem = microgen.phase.BasicGeometry(number=0, shape='polyhedron',
                                        xc=0, yc=0, zc=0,
                                        psi=0, theta=0, phi=0,
                                        param_geom={"dic": dic})
    elem.generate()



    fake = microgen.phase.BasicGeometry(number=0, shape='fake',
                                  xc=0.5, yc=0.5, zc=0.5,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"geom": 0})
    pytest.raises(ValueError, fake.generate)

    microgen.phase.BasicGeometry.__cmp__(elem, fake)
    
    raster = microgen.rasterShapeList(cqShapeList=[ellipsoid], rve=rve, grid=[5, 5, 5])

    compound = cq.Compound.makeCompound(raster[0])
    cq.exporters.export(compound, 'tests/data/compound.step')

    microgen.mesh(mesh_file='tests/data/compound.step', listPhases=raster[1], size=0.03, order=1, output_file="tests/data/compound.msh")

def test_tpms():
    rve = microgen.Rve(dim_x=1, dim_y=1, dim_z=1, size_mesh=0.03)

    elem = microgen.BasicGeometry(number=0, shape='tpms',
                                  xc=0.5, yc=0.5, zc=0.5,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"surface_function": microgen.shape.tpms.gyroid,
                                              "type_part": "skeletal",
                                              "thickness": 0.1},
                                  path_data='tests/data/tpms')
    elem.geometry.createSurfaces(rve=rve,
                                 sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                                 path_data='tests/data/tpms')
    elem.generate(rve=rve)

    elem = microgen.BasicGeometry(number=0, shape='tpms',
                                  xc=0.5, yc=0.5, zc=0.5,
                                  psi=0, theta=0, phi=0,
                                  param_geom={"surface_function": microgen.shape.tpms.schwarzD,
                                              "type_part": "sheet",
                                              "thickness": 0.3},
                                  path_data='tests/data/tpms')
    elem.geometry.createSurfaces(rve=rve,
                                 sizeMesh=0.03, minFacetAngle=20., maxRadius=0.03,
                                 path_data='tests/data/tpms')
    elem.generate(rve=rve)
    
    with pytest.raises(ValueError):
        elem.geometry.createTpms(path_data="tests/data/tpms", rve=None)
    with pytest.raises(ValueError):
        elem.geometry.createTpms(path_data="fake", rve=rve)

    bounding_sphere_radius = 1.1*np.sqrt((0.5*rve.dx)**2 + (0.5*rve.dy)**2 + (0.5*rve.dz)**2)
    generator = microgen.shape.Generator(rve=rve, height=0.3, surface_function=microgen.shape.tpms.gyroid)
    assert generator.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    assert abs(generator.eval([1, 1, 1]) - 0.3) < 1.e-5

    height = 0.2
    assert microgen.shape.tpms.schwarzP(0, 0, 0, height) == 3 + height
    assert microgen.shape.tpms.schwarzD(0, 0, 0, height) == 0 + 0 + 0 + 0 + height
    assert microgen.shape.tpms.neovius(0, 0, 0, height) == (3 + 1 + 1) + (4*1*1*1) + height
    assert microgen.shape.tpms.schoenIWP(0, 0, 0, height) == 2*(1+1+1) - (1+1+1) + height
    assert microgen.shape.tpms.schoenFRD(0, 0, 0, height) == 4 - (1+1+1) + height
    assert microgen.shape.tpms.fischerKochS(0, 0, 0, height) == 0 + 0 + 0 + height
    assert microgen.shape.tpms.pmy(0, 0, 0, height) == 2 + 0 + 0 + 0 + height
    assert microgen.shape.tpms.honeycomb(0, 0, 0, height) == 0 + 0 + 1 + height
    assert microgen.shape.tpms.gyroid(0, 0, 0, height) == 1


if __name__ == "__main__":
    test_shapes()
    test_tpms()