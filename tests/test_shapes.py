import os
import pytest
from inspect import getmembers, isfunction

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


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal", "sheet"])
@pytest.mark.parametrize("surface", [func[0] for func in getmembers(microgen.surface_functions, isfunction)])
def test_tpms_given_cadquery_vtk_shapes_volume_must_be_equivalent(surface: str, type_part: str):
    # Arrange
    tpms = microgen.Tpms(
        surface_function=getattr(microgen.surface_functions, surface),
        offset=1.0,
        resolution=30,
    )

    # Act
    shape = tpms.generate(type_part=type_part, smoothing=0, verbose=False)
    shape_check = tpms.generateVtk(type_part=type_part)

    # Assert
    assert np.isclose(shape.Volume(), np.abs(shape_check.volume), rtol=1e-2)


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal"])
def test_tpms_given_cadquery_vtk_zero_offset_skeletals_volume_must_be_equivalent(type_part: str):
    # Arrange
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=0.0,
    )

    # Act
    shape = tpms.generate(type_part=type_part, smoothing=0, verbose=False)
    shape_check = tpms.generateVtk(type_part=type_part)

    # Assert
    assert np.isclose(shape.Volume(), np.abs(shape_check.volume), rtol=1e-2)

def test_tpms_given_sum_volume_must_be_cube_volume():
    # Arrange
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=1.0,
        repeat_cell=2,
        cell_size=3.0,
    )

    # Act
    volume = np.abs(tpms.sheet.volume + tpms.lower_skeletal.volume + tpms.upper_skeletal.volume)
    cube_volume = np.prod(tpms.repeat_cell) * np.prod(tpms.cell_size)

    # Assert
    assert np.isclose(volume, cube_volume, rtol=1e-2)

@pytest.mark.parametrize("coord_sys_tpms", [microgen.CylindricalTpms, microgen.SphericalTpms])
def test_tpms_given_coord_system_tpms_volumes_must_be_greater_than_zero_and_lower_than_grid_volume(coord_sys_tpms):
    tpms = coord_sys_tpms(
        radius=1.0,
        surface_function=microgen.surface_functions.gyroid,
        offset=1.0,
    )

    assert 0 < tpms.sheet.extract_surface().volume < np.abs(tpms.grid.volume)
    assert 0 < tpms.lower_skeletal.extract_surface().volume < np.abs(tpms.grid.volume)
    assert 0 < tpms.upper_skeletal.extract_surface().volume < np.abs(tpms.grid.volume)

def test_tpms_given_variable_thickness_tpms_volumes_must_be_greater_than_zero_and_lower_than_grid_volume():
    def graded_density(x, y, z):
        return x

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=graded_density,
    )

    assert 0 < tpms.sheet.extract_surface().volume < np.abs(tpms.grid.volume)
    assert 0 < tpms.lower_skeletal.extract_surface().volume < np.abs(tpms.grid.volume)
    assert 0 < tpms.upper_skeletal.extract_surface().volume < np.abs(tpms.grid.volume)

def test_tpms_given_incorrect_parameters_must_raise_errors():
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=0.5,
    )
    with pytest.raises(ValueError):
        tpms.generateVtk(type_part="fake")
    with pytest.raises(ValueError):
        tpms.generate(type_part="fake")

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