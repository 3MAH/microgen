from microgen import periodic, Rve, shape, Phase

import pytest


def generate_sphere(x, y, z, rve):
    elem = shape.sphere.Sphere(center=(x, y, z), radius=0.1)
    phase = Phase(shape=elem.generate())
    periodic(phase=phase, rve=rve)


@pytest.mark.filterwarnings("ignore:Object intersecting")
def test_periodic():
    rve = Rve(dim_x=1, dim_y=1, dim_z=1, center=(0.5, 0.5, 0.5))

    # test x- and x+ faces intersected
    elem = shape.capsule.Capsule(center=(0.5, 0, 0.5), height=1, radius=0.1)
    phase = Phase(shape=elem.generate())
    periodic(phase=phase, rve=rve)

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


if __name__ == "__main__":
    test_periodic()
