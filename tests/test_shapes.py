import microgen

import cadquery as cq
import numpy as np

import os
import pytest


def test_shapes():
    os.makedirs("tests/data", exist_ok=True)  # if data folder doesn't exist yet

    rve = microgen.rve.Rve(dim_x=1, dim_y=1, dim_z=1)

    elem = microgen.shape.newGeometry(
        shape="Ellipsoid", param_geom={"a_x": 0.15, "a_y": 0.31, "a_z": 0.4}
    )
    elem = microgen.shape.ellipsoid.Ellipsoid(a_x=0.15, a_y=0.31, a_z=0.4)
    ellipsoid = elem.generate()
    elem.generateVtk()
    phase = microgen.phase.Phase(shape=ellipsoid)
    phase.centerOfMass
    phase.centerOfMass
    phase.getCenterOfMass(compute=False)
    phase.inertiaMatrix
    phase.inertiaMatrix
    phase.getInertiaMatrix(compute=False)
    phase.solids
    phase.shape
    phase.translate((1, 0, 0))
    phase.translate(np.array([0, 1, 1]))
    phase.rescale(1.5)
    phase.repeat(rve, (1, 2, 1))
    phase.rasterize(rve, [2, 2, 2], phasePerRaster=True)
    phase.rasterize(rve, [2, 2, 2], phasePerRaster=False)

    void_phase = microgen.phase.Phase()
    void_phase.shape
    void_phase.solids

    elem = microgen.shape.newGeometry(shape="Sphere", param_geom={"radius": 0.15})
    elem = microgen.shape.sphere.Sphere(radius=0.15)
    elem.generate()
    elem.generateVtk()

    elem = microgen.shape.newGeometry(
        shape="Box", param_geom={"dim_x": 0.15, "dim_y": 0.31, "dim_z": 0.4}
    )
    elem = microgen.shape.box.Box(dim_x=0.15, dim_y=0.31, dim_z=0.4)
    elem.generate()
    elem.generateVtk()

    elem = microgen.shape.newGeometry(
        shape="Capsule", param_geom={"height": 0.5, "radius": 0.1}
    )
    elem = microgen.shape.capsule.Capsule(height=0.5, radius=0.1)
    elem.generate()
    elem.generateVtk()

    elem = microgen.shape.newGeometry(
        shape="Cylinder", param_geom={"height": 0.5, "radius": 0.1}
    )
    elem = microgen.shape.cylinder.Cylinder(height=0.5, radius=0.1)
    elem.generate()
    elem.generateVtk()

    elem = microgen.shape.newGeometry(
        shape="ExtrudedPolygon",
        param_geom={"listCorners": [(0, 0), (0, 1), (1, 1), (1, 0)], "height": 0.3},
    )
    elem = microgen.shape.extrudedPolygon.ExtrudedPolygon(
        listCorners=[(0, 0), (0, 1), (1, 1), (1, 0)], height=0.3
    )
    elem.generate()
    elem.generateVtk()

    elem = microgen.shape.polyhedron.Polyhedron()  # default shape = tetrahedron
    elem.generate()
    elem.generateVtk()
    dic = microgen.shape.polyhedron.read_obj(
        "examples/BasicShapes/platon/tetrahedron.obj"
    )
    microgen.shape.newGeometry(shape="Polyhedron", param_geom={"dic": dic})

    with pytest.raises(ValueError):
        microgen.shape.newGeometry(shape="fake", param_geom={"fake": 0})

    raster = microgen.operations.rasterPhase(
        phase=phase, rve=rve, grid=[5, 5, 5], phasePerRaster=False
    )

    raster = microgen.operations.rasterPhase(
        phase=phase, rve=rve, grid=[5, 5, 5]
    )

    compound = cq.Compound.makeCompound([solid for phase in raster for solid in phase.solids])
    cq.exporters.export(compound, "tests/data/compound.step")

    microgen.mesh(
        mesh_file="tests/data/compound.step",
        listPhases=raster,
        size=0.03,
        order=1,
        output_file="tests/data/compound.msh",
    )


def test_tpms():
    elem = microgen.shape.newGeometry(
        shape="tpms",
        param_geom={
            "surface_function": microgen.shape.tpms.gyroid,
            "type_part": "skeletal",
            "thickness": 0.075,
            "cell_size": 1,
            "repeat_cell": 1,
        },
    )
    elem.generate()
    elem.generateVtk()

    elem = microgen.shape.tpms.Tpms(
        center=(0.5, 0.5, 0.5),
        surface_function=microgen.shape.tpms.schwarzD,
        type_part="sheet",
        thickness=0.05,
        cell_size=(1, 2, 1),
        repeat_cell=(2, 1, 1),
    )
    elem.generate()
    elem.generateSurface(isovalue=0.1)
    elem.generateSurfaceVtk()

    with pytest.raises(ValueError):
        microgen.shape.tpms.Tpms(
            center=(0.5, 0.5, 0.5),
            surface_function=microgen.shape.tpms.schwarzD,
            type_part="fake",
            thickness=0.3,
        )

    assert microgen.shape.tpms.schwarzP(0, 0, 0) == 3
    assert microgen.shape.tpms.schwarzD(0, 0, 0) == 0 + 0 + 0 + 0
    assert (
        microgen.shape.tpms.neovius(0, 0, 0)
        == (3 + 1 + 1) + (4 * 1 * 1 * 1)
    )
    assert (
        microgen.shape.tpms.schoenIWP(0, 0, 0)
        == 2 * (1 + 1 + 1) - (1 + 1 + 1)
    )
    assert microgen.shape.tpms.schoenFRD(0, 0, 0) == 4 - (1 + 1 + 1)
    assert microgen.shape.tpms.fischerKochS(0, 0, 0) == 0 + 0 + 0
    assert microgen.shape.tpms.pmy(0, 0, 0) == 2 + 0 + 0 + 0
    assert microgen.shape.tpms.honeycomb(0, 0, 0) == 0 + 0 + 1
    assert microgen.shape.tpms.gyroid(0, 0, 0) == 0


if __name__ == "__main__":
    test_shapes()
    test_tpms()
