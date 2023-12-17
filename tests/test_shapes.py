import os
from inspect import getmembers, isfunction
from typing import Literal, Type, Union

import cadquery as cq
import numpy as np
import pytest

import microgen
from microgen.shape import (
    Box,
    Capsule,
    Cylinder,
    Ellipsoid,
    ExtrudedPolygon,
    Polyhedron,
    Sphere,
)


def test_shapes():
    rve = microgen.Rve(dim=1)

    elem = microgen.newGeometry(
        shape="Ellipsoid", param_geom={"a_x": 0.15, "a_y": 0.31, "a_z": 0.4}
    )
    ellipsoid = elem.generate()
    phase = microgen.Phase(shape=ellipsoid)
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

    void_phase = microgen.Phase()
    void_phase.shape
    void_phase.solids

    elem = microgen.newGeometry(shape="Sphere", param_geom={"radius": 0.15})
    elem = microgen.newGeometry(
        shape="Box", param_geom={"dim_x": 0.15, "dim_y": 0.31, "dim_z": 0.4}
    )
    elem = microgen.newGeometry(
        shape="Capsule", param_geom={"height": 0.5, "radius": 0.1}
    )
    elem = microgen.newGeometry(
        shape="Cylinder", param_geom={"height": 0.5, "radius": 0.1}
    )
    elem = microgen.newGeometry(
        shape="ExtrudedPolygon",
        param_geom={"listCorners": [(0, 0), (0, 1), (1, 1), (1, 0)], "height": 0.3},
    )
    dic = microgen.shape.polyhedron.read_obj(
        "examples/BasicShapes/platon/tetrahedron.obj"
    )
    microgen.newGeometry(shape="Polyhedron", param_geom={"dic": dic})

    with pytest.raises(ValueError):
        microgen.newGeometry(shape="fake", param_geom={"fake": 0})


@pytest.mark.parametrize(
    "shape",
    [Box, Capsule, Cylinder, Ellipsoid, ExtrudedPolygon, Polyhedron, Sphere],
)
def test_shape_cad_and_vtk_volume_must_correspond(shape: Type[microgen.BasicGeometry]):
    geom = shape()
    shape_cad = geom.generate()
    shape_vtk = geom.generateVtk()

    volume_cad = shape_cad.Volume()

    assert volume_cad > 0
    assert np.isclose(volume_cad, shape_vtk.volume, rtol=1e-2)


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal", "sheet"])
def test_tpms_given_cadquery_vtk_shapes_volume_must_be_equivalent(
    type_part: Literal["sheet", "lower skeletal", "upper skeletal"],
):
    # Arrange
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        density=0.3,
    )

    # Act
    shape_cadquery = tpms.generate(type_part=type_part, smoothing=0, verbose=False)
    shape_vtk = tpms.generateVtk(type_part=type_part)

    # Assert
    assert np.isclose(shape_cadquery.Volume(), np.abs(shape_vtk.volume), rtol=1e-2)


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal"])
def test_tpms_given_cadquery_vtk_zero_offset_skeletals_volume_must_be_equivalent(
    type_part: Literal["lower skeletal", "upper skeletal"],
):
    # Arrange
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.schwarzP,
        offset=0.0,
    )

    # Act
    shape_cadquery = tpms.generate(type_part=type_part, smoothing=0, verbose=False)
    shape_vtk = tpms.generateVtk(type_part=type_part)

    # Assert
    assert np.isclose(shape_cadquery.Volume(), np.abs(shape_vtk.volume), rtol=1e-2)


@pytest.mark.parametrize(
    "surface", [func[0] for func in getmembers(microgen.surface_functions, isfunction)]
)
@pytest.mark.parametrize("repeat_cell", [2, (2, 1, 3)])
@pytest.mark.parametrize("cell_size", [3.0, (0.5, 1.5, 1.0)])
def test_tpms_given_sum_volume_must_be_cube_volume(
    surface: str,
    repeat_cell: Union[int, tuple[int, int, int]],
    cell_size: Union[float, tuple[float, float, float]],
):
    # Arrange
    tpms = microgen.Tpms(
        surface_function=getattr(microgen.surface_functions, surface),
        offset=1.0,
        repeat_cell=repeat_cell,
        cell_size=cell_size,
    )

    # Act
    volume = np.abs(
        tpms.sheet.volume + tpms.lower_skeletal.volume + tpms.upper_skeletal.volume
    )
    cube_volume = np.prod(tpms.repeat_cell) * np.prod(tpms.cell_size)

    # Assert
    assert np.isclose(volume, cube_volume, rtol=1e-2)


@pytest.mark.parametrize(
    "surface", [func[0] for func in getmembers(microgen.surface_functions, isfunction)]
)
@pytest.mark.parametrize("density", [0.3, 0.7])
def test_tpms_given_density_must_match_computed_density(
    surface: str,
    density: float,
):
    # Arrange
    tpms = microgen.Tpms(
        surface_function=getattr(microgen.surface_functions, surface), density=density
    )

    # Act
    computed_density = tpms.generateVtk(type_part="sheet").volume / tpms.grid.volume

    # Assert
    assert np.isclose(computed_density, density, rtol=1e-2)


@pytest.mark.parametrize(
    "coord_sys_tpms", [microgen.CylindricalTpms, microgen.SphericalTpms]
)
def test_tpms_given_coord_system_tpms_volumes_must_be_greater_than_zero_and_lower_than_grid_volume(
    coord_sys_tpms: Type[microgen.Tpms],
):
    tpms = coord_sys_tpms(
        radius=1.0,
        surface_function=microgen.surface_functions.gyroid,
        density=0.2,
    )

    assert 0 < tpms.sheet.extract_surface().volume < np.abs(tpms.grid.volume)
    assert 0 < tpms.lower_skeletal.extract_surface().volume < np.abs(tpms.grid.volume)
    assert 0 < tpms.upper_skeletal.extract_surface().volume < np.abs(tpms.grid.volume)


@pytest.mark.parametrize(
    "coord_sys_tpms,repeat_cell_zero,repeat_cell_max",
    [
        (microgen.CylindricalTpms, (1, 0, 1), (1, 100, 1)),
        (microgen.SphericalTpms, (1, 0, 0), (1, 100, 100)),
    ],
)
def test_tpms_given_zero_and_max_repeat_cell_values_volumes_must_correspond_and_be_between_zero_and_grid_volume(
    coord_sys_tpms: Type[microgen.Tpms],
    repeat_cell_zero: tuple[int, int, int],
    repeat_cell_max: tuple[int, int, int],
):
    tpms_repeat_zero = coord_sys_tpms(
        radius=1.0,
        surface_function=microgen.surface_functions.gyroid,
        offset=1.0,
        repeat_cell=repeat_cell_zero,
    )

    tpms_repeat_max = coord_sys_tpms(
        radius=1.0,
        surface_function=microgen.surface_functions.gyroid,
        offset=1.0,
        repeat_cell=repeat_cell_max,
    )

    assert np.abs(tpms_repeat_zero.grid.volume) == np.abs(tpms_repeat_max.grid.volume)

    assert (
        0
        < tpms_repeat_zero.sheet.extract_surface().volume
        == tpms_repeat_max.sheet.extract_surface().volume
        < np.abs(tpms_repeat_zero.grid.volume)
    )
    assert (
        0
        < tpms_repeat_zero.lower_skeletal.extract_surface().volume
        == tpms_repeat_max.lower_skeletal.extract_surface().volume
        < np.abs(tpms_repeat_zero.grid.volume)
    )
    assert (
        0
        < tpms_repeat_zero.upper_skeletal.extract_surface().volume
        == tpms_repeat_max.upper_skeletal.extract_surface().volume
        < np.abs(tpms_repeat_zero.grid.volume)
    )


def test_tpms_given_generate_surface_must_not_be_empty():
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=1.0,
        density=0.2,
    )

    surface = tpms.generate(type_part="surface")
    assert np.any(surface.Vertices())
    assert np.any(surface.Faces())
    assert not surface.Closed()


def test_tpms_given_variable_offset_cadquery_and_vtk_volumes_must_correspond():
    def variable_offset(x: np.ndarray, y: np.ndarray, z: np.ndarray):
        return x + 1.5

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=variable_offset,
    )

    shape_cadquery = tpms.generate(type_part="sheet", smoothing=0, verbose=True)
    shape_vtk = tpms.generateVtk(type_part="sheet")

    assert np.isclose(shape_cadquery.Volume(), np.abs(shape_vtk.volume), rtol=1e-2)


@pytest.mark.parametrize("param", [0.3, 4.0])
def test_tpms_given_variable_offset_out_of_limits_with_cadquery_must_raise_ValueError(
    param: float,
):
    def variable_offset(x: np.ndarray, y: np.ndarray, z: np.ndarray):
        return x + param  # x ∈ [-0.5, 0.5]

    # offset must be in [0, 2 * max(gyroid)]
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=variable_offset,
    )

    with pytest.raises((ValueError, NotImplementedError)):
        tpms.generate(type_part="sheet")


def test_tpms_generate_given_wrong_type_part_parameter_must_raise_ValueError():
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=0.5,
    )
    with pytest.raises(ValueError):
        tpms.generateVtk(type_part="fake")
    with pytest.raises(ValueError):
        tpms.generate(type_part="fake")


def test_tpms_given_wrong_cell_size_parameter_must_raise_ValueError():
    with pytest.raises(ValueError):
        microgen.Tpms(
            surface_function=microgen.surface_functions.gyroid,
            cell_size=(1, 1),
        )


def test_tpms_given_wrong_repeat_cell_parameter_must_raise_ValueError():
    with pytest.raises(ValueError):
        microgen.Tpms(
            surface_function=microgen.surface_functions.gyroid,
            repeat_cell=(1, 1, 1, 1),
        )


def test_tpms_given_wrong_density_parameter_must_raise_ValueError():
    with pytest.raises(ValueError):
        microgen.Tpms(
            surface_function=microgen.surface_functions.gyroid,
            density=0.0,
        )

    with pytest.raises(ValueError):
        microgen.Tpms.offset_from_density(
            surface_function=microgen.surface_functions.gyroid,
            density="fake",
            part_type="sheet",
        )

    with pytest.raises(ValueError):
        microgen.Tpms.offset_from_density(
            surface_function=microgen.surface_functions.gyroid,
            density=1.0,
            part_type="lower skeletal",
        )


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal", "sheet"])
def test_tpms_given_density_must_generate_tpms_with_correct_volume(type_part: str):
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        density=0.2,
    )

    part = tpms.generateVtk(type_part=type_part)
    assert np.isclose(part.volume, tpms.grid.volume * 0.2, rtol=1e-2)


@pytest.mark.parametrize("part_type", ["lower skeletal", "upper skeletal", "sheet"])
def test_tpms_given_max_density_must_return_corresponding_offset(
    part_type: Literal["sheet", "lower skeletal", "upper skeletal"]
):
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
    )

    max_offset = microgen.Tpms.offset_from_density(
        surface_function=microgen.surface_functions.gyroid,
        density="max",
        part_type=part_type,
    )
    expected_max_offset = (
        2.0 * np.max(tpms.grid["surface"]) if part_type == "sheet" else 0.0
    )
    assert max_offset == expected_max_offset


def test_tpms_given_100_percent_sheet_density_must_return_a_cube():
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        density=1.0,
    )

    assert np.isclose(tpms.sheet.volume, tpms.grid.volume, rtol=1.0e-9)


def test_tpms_given_density_must_return_corresponding_offset():
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
    )

    offset = microgen.Tpms.offset_from_density(
        surface_function=microgen.surface_functions.gyroid,
        density=0.5,
        part_type="sheet",
    )
    assert 0 < offset < 2.0 * np.max(tpms.grid["surface"])


def test_tpms_given_property_must_return_the_same_value():
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        density=0.2,
    )
    skeletals = tpms.skeletals

    assert tpms.upper_skeletal == skeletals[0]
    assert tpms.lower_skeletal == skeletals[1]
    assert tpms.sheet == tpms.sheet
    assert tpms.generateVtk(type_part="surface") == tpms.surface


def test_tpms_given_surface_must_not_be_empty():
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
    )

    assert np.any(tpms.surface.points)
    assert np.any(tpms.surface.faces)


def test_tpms_given_negative_offset_for_skeletal_must_work_with_vtk_and_raise_error_with_cadquery():
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=-1.0,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate(type_part="lower skeletal")

    sheet = tpms.generateVtk(type_part="lower skeletal").extract_surface()
    assert 0.0 < sheet.volume < np.abs(tpms.grid.volume)

    def including_negative_values(x: np.ndarray, y: np.ndarray, z: np.ndarray):
        return x

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=including_negative_values,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate(type_part="lower skeletal")

    sheet = tpms.generateVtk(type_part="lower skeletal").extract_surface()
    assert 0.0 < sheet.volume < np.abs(tpms.grid.volume)


def test_tpms_given_negative_offset_for_sheet_must_work_with_vtk_and_raise_error_with_cadquery():
    def all_negative(x: np.ndarray, y: np.ndarray, z: np.ndarray):
        return -1.0 + x

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=all_negative,
    )
    with pytest.raises(ValueError):
        tpms.generate(type_part="sheet")

    assert tpms.generateVtk(type_part="sheet").volume == 0.0

    def including_negative_values(x: np.ndarray, y: np.ndarray, z: np.ndarray):
        return x

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=including_negative_values,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate(type_part="sheet")

    sheet = tpms.generateVtk(type_part="sheet").extract_surface()
    assert 0.0 < sheet.volume < np.abs(tpms.grid.volume)
