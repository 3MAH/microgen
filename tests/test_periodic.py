import microgen

def generate_sphere(x, y, z, rve):
    elem = microgen.BasicGeometry(shape='sphere',
                                  xc=x, yc=y, zc=z,
                                  param_geom={"radius": 0.1})
    phase = elem.generate()
    microgen.periodic(cqshape=phase, rve=rve)


def test_periodic():
    rve = microgen.Rve(dim_x=1, dim_y=1, dim_z=1, size_mesh=0.03)

    # test x- and x+ faces intersected 
    elem = microgen.BasicGeometry(shape='capsule',
                                  xc=0.5, yc=0, zc=0.5,
                                  param_geom={"height": 1,
                                              "radius": 0.1})
    phase = elem.generate()
    microgen.periodic(cqshape=phase, rve=rve)

    # test no intersection
    generate_sphere(x=0.5, y=0.5, z=0.5, rve=rve)

    # face
    generate_sphere(x=0, y=0.5, z=0.5, rve=rve)
    generate_sphere(x=1, y=0.5, z=0.5, rve=rve)
    generate_sphere(x=0.5, y=0, z=0.5, rve=rve)
    generate_sphere(x=0.5, y=1, z=0.5, rve=rve)
    generate_sphere(x=0.5, y=0.5, z=0, rve=rve)
    generate_sphere(x=0.5, y=0.5, z=1, rve=rve)

    # edge
    generate_sphere(x=0, y=0, z=0.5, rve=rve)
    generate_sphere(x=0, y=1, z=0.5, rve=rve)
    generate_sphere(x=0, y=0.5, z=0, rve=rve)
    generate_sphere(x=0, y=0.5, z=1, rve=rve)
    generate_sphere(x=1, y=0, z=0.5, rve=rve)
    generate_sphere(x=1, y=1, z=0.5, rve=rve)
    generate_sphere(x=1, y=0.5, z=0, rve=rve)
    generate_sphere(x=1, y=0.5, z=1, rve=rve)
    generate_sphere(x=0.5, y=0, z=0, rve=rve)
    generate_sphere(x=0.5, y=0, z=1, rve=rve)
    generate_sphere(x=0.5, y=1, z=0, rve=rve)
    generate_sphere(x=0.5, y=1, z=1, rve=rve)

    # corner
    generate_sphere(x=0, y=0, z=0, rve=rve)
    generate_sphere(x=0, y=0, z=1, rve=rve)
    generate_sphere(x=0, y=1, z=0, rve=rve)
    generate_sphere(x=0, y=1, z=1, rve=rve)
    generate_sphere(x=1, y=0, z=0, rve=rve)
    generate_sphere(x=1, y=0, z=1, rve=rve)
    generate_sphere(x=1, y=1, z=0, rve=rve)
    generate_sphere(x=1, y=1, z=1, rve=rve)



