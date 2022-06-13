import microgen


def test_misc():
    microgen.parseNeper("examples/3Doperations/Voronoi/test1")

    microgen.removeEmptyLines("fake_file.txt")

    microgen.MatSection(
        number=0,
        name="test",
        umat_name="test",
        psi_mat=0,
        theta_mat=0,
        phi_mat=0,
        nprops=0,
        nstatev=0,
        props=[],
    )


def test_operations():

    elem = microgen.shape.Box(center=(0.5, 0.5, 0.5), dim_x=1, dim_y=1, dim_z=1)
    shape1 = elem.generate()
    phase1 = microgen.Phase(shape=shape1)

    elem = microgen.shape.Box(center=(0, 0, 0), dim_x=0.5, dim_y=0.5, dim_z=0.5)
    shape2 = elem.generate()

    microgen.cutPhaseByShapeList(phaseToCut=phase1, cqShapeList=[shape2])

    microgen.cutPhasesByShape(phaseList=[phase1], cut_obj=shape2)

    rve = microgen.Rve(dim_x=1, dim_y=1, dim_z=1)
    microgen.repeatGeometry(phase1, rve, grid=[2, 2, 2])

    #  EXTERNAL
    #  runNeper

    #  MATERIAL
    #  readSections


if __name__ == "__main__":
    test_misc()
    test_operations()
