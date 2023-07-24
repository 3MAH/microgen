import os
import pytest

import cadquery as cq
import numpy as np

import microgen


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
            "surface_function": microgen.shape.surface_functions.gyroid,
            "offset": 0.3,
            "cell_size": (1, 2, 1),
            "repeat_cell": (1, 2, 1),
            "resolution": 10,
        },
    )
    _ = elem.surface
    _ = elem.surface
    _ = elem.skeletals
    sheet = elem.sheet
    _ = elem.sheet
    assert 0 < sheet.volume < np.prod(elem.cell_size) * np.prod(elem.repeat_cell)
    _ = elem.generateVtk(type_part="lower skeletal")
    _ = elem.generateVtk(type_part="upper skeletal")
    _ = elem.generateVtk(type_part="surface")
    _ = elem.generate(type_part="sheet")

    def graded_density(x, y, z):
        return x
    elem = microgen.shape.tpms.CylindricalTpms(
        radius=1.0,
        surface_function=microgen.shape.surface_functions.schwarzD,
        offset=graded_density,
        phase_shift=(0.1, 0.2, 0.3),
        cell_size=(1, 2, 1),
        repeat_cell=(2, 1, 1),
        resolution=20,
    )
    # elem.generate(type_part="surface")
    _ = elem.upper_skeletal
    _ = elem.upper_skeletal
    _ = elem.generateVtk(type_part="sheet")

    elem = microgen.shape.tpms.CylindricalTpms(
        radius=1.0,
        surface_function=microgen.shape.surface_functions.schwarzD,
        repeat_cell=(1, 10, 1),
    )
    _ = elem.lower_skeletal
    _ = elem.lower_skeletal

    elem = microgen.shape.tpms.SphericalTpms(
        radius=1.0,
        surface_function=microgen.shape.surface_functions.schwarzD,
        repeat_cell=(1, 10, 10),
    )
    _ = elem.surface
    _ = elem.surface

    with pytest.raises(ValueError):
        _ = elem.generateVtk(type_part="fake")

    with pytest.raises(ValueError):
        microgen.shape.tpms.Tpms(
            surface_function=microgen.shape.surface_functions.gyroid,
            cell_size=(1, 1),
        )

    with pytest.raises(ValueError):
        microgen.shape.tpms.Tpms(
            surface_function=microgen.shape.surface_functions.gyroid,
            repeat_cell=(1, 1, 1, 1),
        )

    assert microgen.shape.surface_functions.schwarzP(0, 0, 0) == 3
    assert microgen.shape.surface_functions.schwarzD(0, 0, 0) == 0 + 0 + 0 + 0
    assert (
        microgen.shape.surface_functions.neovius(0, 0, 0)
        == (3 + 1 + 1) + (4 * 1 * 1 * 1)
    )
    assert (
        microgen.shape.surface_functions.schoenIWP(0, 0, 0)
        == 2 * (1 + 1 + 1) - (1 + 1 + 1)
    )
    assert microgen.shape.surface_functions.schoenFRD(0, 0, 0) == 4 - (1 + 1 + 1)
    assert microgen.shape.surface_functions.fischerKochS(0, 0, 0) == 0 + 0 + 0
    assert microgen.shape.surface_functions.pmy(0, 0, 0) == 2 + 0 + 0 + 0
    assert microgen.shape.surface_functions.honeycomb(0, 0, 0) == 0 + 0 + 1
    assert microgen.shape.surface_functions.gyroid(0, 0, 0) == 0
    assert round(microgen.shape.surface_functions.lidinoid(0, 0, 0), 2) == -1.2
    assert round(microgen.shape.surface_functions.split_p(0, 0, 0), 2) == -1.8


if __name__ == "__main__":
    test_shapes()
    test_tpms()
