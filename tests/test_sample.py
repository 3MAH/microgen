import microgen

from sys import platform

def test_misc():
    if platform != 'win32':
        microgen.Neper.run(filename='tests/data/neper.tess', nbCell=2, dimCube=(1, 1, 1))
        microgen.parseNeper("tests/data/neper")
        microgen.Neper.generateVoronoiFromTessFile('tests/data/neper.tess')


def test_operations():
    elem = microgen.shape.Box(center=(0.5, 0.5, 0.5), dim_x=1, dim_y=1, dim_z=1)
    shape1 = elem.generate()
    phase1 = microgen.Phase(shape=shape1)

    elem = microgen.shape.Box(center=(0, 0, 0), dim_x=0.5, dim_y=0.5, dim_z=0.5)
    shape2 = elem.generate()
    microgen.rescale(shape2, 2.)

    microgen.cutPhaseByShapeList(phaseToCut=phase1, cqShapeList=[shape2])

    microgen.cutPhasesByShape(phaseList=[phase1], cut_obj=shape2)

    rve = microgen.Rve(dim_x=1, dim_y=1, dim_z=1, center=(0.5, 0.5, 0.5))
    microgen.repeatShape(shape1, rve, grid=[2, 2, 2])


if __name__ == "__main__":
    test_misc()
    test_operations()
