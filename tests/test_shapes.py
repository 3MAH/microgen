import microgen

import cadquery as cq
import numpy as np

import pytest


def test_shapes():
    rve = microgen.rve.Rve(dim_x=1, dim_y=1, dim_z=1)

    elem = microgen.shape.newGeometry(
        shape="Ellipsoid", param_geom={"a_x": 0.15, "a_y": 0.31, "a_z": 0.4}
    )
    elem = microgen.shape.ellipsoid.Ellipsoid(a_x=0.15, a_y=0.31, a_z=0.4)
    ellipsoid = elem.generate()
    phase = microgen.phase.Phase(shape=ellipsoid)
    phase.centerOfMass
    phase.centerOfMass
    phase.inertiaMatrix
    phase.inertiaMatrix
    phase.solids
    phase.shape

    microgen.phase.Phase()

    elem = microgen.shape.newGeometry(shape="Sphere", param_geom={"radius": 0.15})
    elem = microgen.shape.sphere.Sphere(radius=0.15)
    elem.generate()

    elem = microgen.shape.newGeometry(
        shape="Box", param_geom={"dim_x": 0.15, "dim_y": 0.31, "dim_z": 0.4}
    )
    elem = microgen.shape.box.Box(dim_x=0.15, dim_y=0.31, dim_z=0.4)
    elem.generate()

    elem = microgen.shape.newGeometry(
        shape="Capsule", param_geom={"height": 0.5, "radius": 0.1}
    )
    elem = microgen.shape.capsule.Capsule(height=0.5, radius=0.1)
    elem.generate()

    elem = microgen.shape.newGeometry(
        shape="Cylinder", param_geom={"height": 0.5, "radius": 0.1}
    )
    elem = microgen.shape.cylinder.Cylinder(height=0.5, radius=0.1)
    elem.generate()

    elem = microgen.shape.newGeometry(
        shape="ExtrudedPolygon",
        param_geom={"listCorners": [(0, 0), (0, 1), (1, 1), (1, 0)], "height": 0.3},
    )
    elem = microgen.shape.extrudedPolygon.ExtrudedPolygon(
        listCorners=[(0, 0), (0, 1), (1, 1), (1, 0)], height=0.3
    )
    elem.generate()

    # elem = microgen.shape.polyhedron.Polyhedron()  # default shape = tetrahedron
    # elem.generate()
    dic = microgen.shape.polyhedron.read_obj(
        "examples/BasicShapes/platon/tetrahedron.obj"
    )
    microgen.shape.newGeometry(shape="Polyhedron", param_geom={"dic": dic})

    with pytest.raises(ValueError):
        microgen.shape.newGeometry(shape="fake", param_geom={"fake": 0})

    raster = microgen.operations.rasterShapeList(
        cqShapeList=[ellipsoid], rve=rve, grid=[5, 5, 5]
    )

    compound = cq.Compound.makeCompound(raster[0])
    cq.exporters.export(compound, "tests/data/compound.step")

    listPhases = [microgen.Phase(solids=solidList) for solidList in raster[1]]
    microgen.mesh(
        mesh_file="tests/data/compound.step",
        listPhases=listPhases,
        size=0.03,
        order=1,
        output_file="tests/data/compound.msh",
    )


def test_tpms():
    elem = microgen.shape.newGeometry(
        shape="tpms",
        center=(0.5, 0.5, 0.5),
        param_geom={
            "surface_function": microgen.shape.tpms.gyroid,
            "type_part": "skeletal",
            "thickness": 0.1,
            "cell_size": None,
            "repeat_cell": 1,
            "path_data": "tests/data/tpms1",
        },
    )

    elem = microgen.shape.tpms.Tpms(
        center=(0.5, 0.5, 0.5),
        surface_function=microgen.shape.tpms.schwarzD,
        type_part="sheet",
        thickness=0.3,
        cell_size=2,
        repeat_cell=(2, 1, 1),
        path_data="tests/data/tpms2",
    )
    elem.generate()

    elem = microgen.shape.tpms.Tpms(
        center=(0.5, 0.5, 0.5),
        surface_function=microgen.shape.tpms.schwarzD,
        type_part="sheet",
        thickness=0.3,
        cell_size=2,
        repeat_cell=(2, 1, 1),
        path_data="tests/data/tpms2",
    )
    elem.generateSurface()
    elem.generate()

    with pytest.raises(ValueError):
        microgen.shape.tpms.Tpms(
            center=(0.5, 0.5, 0.5),
            surface_function=microgen.shape.tpms.schwarzD,
            type_part="fake",
            thickness=0.3,
        )

    bounding_sphere_radius = 1.1 * np.sqrt(0.5**2 + 0.5**2 + 0.5**2)
    generator = microgen.shape.tpms.Generator(
        height=0.3, surface_function=microgen.shape.tpms.gyroid
    )
    assert generator.get_bounding_sphere_squared_radius() == bounding_sphere_radius
    assert abs(generator.eval([1, 1, 1]) - 0.3) < 1.0e-5

    height = 0.2
    assert microgen.shape.tpms.schwarzP(0, 0, 0, height) == 3 + height
    assert microgen.shape.tpms.schwarzD(0, 0, 0, height) == 0 + 0 + 0 + 0 + height
    assert (
        microgen.shape.tpms.neovius(0, 0, 0, height)
        == (3 + 1 + 1) + (4 * 1 * 1 * 1) + height
    )
    assert (
        microgen.shape.tpms.schoenIWP(0, 0, 0, height)
        == 2 * (1 + 1 + 1) - (1 + 1 + 1) + height
    )
    assert microgen.shape.tpms.schoenFRD(0, 0, 0, height) == 4 - (1 + 1 + 1) + height
    assert microgen.shape.tpms.fischerKochS(0, 0, 0, height) == 0 + 0 + 0 + height
    assert microgen.shape.tpms.pmy(0, 0, 0, height) == 2 + 0 + 0 + 0 + height
    assert microgen.shape.tpms.honeycomb(0, 0, 0, height) == 0 + 0 + 1 + height
    assert microgen.shape.tpms.gyroid(0, 0, 0, height) == 1


if __name__ == "__main__":
    test_shapes()
    test_tpms()
