"""Tests for the shapes module."""

from __future__ import annotations

from inspect import getmembers, isfunction
from typing import Any, Dict, Literal, Tuple, Type

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


@pytest.mark.parametrize(
    ("shape", "params"),
    [
        (Box, {"dim": (0.15, 0.31, 0.4)}),
        (Capsule, {"height": 0.5, "radius": 0.1}),
        (Cylinder, {"height": 0.5, "radius": 0.1}),
        (Ellipsoid, {"radii": (0.15, 0.31, 0.4)}),
        (
            ExtrudedPolygon,
            {"listCorners": [(0, 0), (0, 1), (1, 1), (1, 0)], "height": 0.3},
        ),
        (
            Polyhedron,
            {
                "dic": microgen.shape.polyhedron.read_obj(
                    "examples/BasicShapes/platon/tetrahedron.obj"
                )
            },
        ),
        (Sphere, {"radius": 0.15}),
    ],
)
def test_new_geometry_and_native_shape_must_have_same_volume(
    shape: Type[microgen.BasicGeometry], params: Dict[str, Any]
):
    """Test if the newGeometry generates the same volume as the native shape."""
    key = shape.__name__
    new_geom = microgen.newGeometry(shape=key, param_geom=params)
    native_shape = shape(**params)

    assert new_geom.generateVtk().volume == native_shape.generateVtk().volume


def test_new_geometry_given_wrong_shape_must_raise_ValueError():
    """Test if the newGeometry raises an error when the shape is not recognized."""
    with pytest.raises(ValueError):
        microgen.newGeometry(shape="fake", param_geom={"fake": 0})


@pytest.mark.parametrize(
    "shape",
    [Box, Capsule, Cylinder, Ellipsoid, ExtrudedPolygon, Polyhedron, Sphere],
)
def test_shape_cad_and_vtk_volume_must_correspond(shape: Type[microgen.BasicGeometry]):
    """Test if the shape generates the same volume with CadQuery and VTK."""
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
    """Test if the TPMS generates the same volume with CadQuery and VTK."""
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
    """Test if the TPMS zero offset generates the same volume with CadQuery and VTK."""
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
    repeat_cell: int | Tuple[int, int, int],
    cell_size: float | Tuple[float, float, float],
):
    """Test if the TPMS sum volume is equal to the cube volume."""
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
    """Test if the TPMS density matches the computed density."""
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
    """Test if the TPMS volumes are between zero and the grid volume."""
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
def test_tpms_given_zero_and_max_repeat_cell_values_volumes_must_correspond(
    coord_sys_tpms: Type[microgen.Tpms],
    repeat_cell_zero: Tuple[int, int, int],
    repeat_cell_max: Tuple[int, int, int],
):
    """Test if the TPMS volumes correspond and are between zero and the grid volume."""
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
    """Test if the TPMS generate surface is not empty."""
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
    """Test if the TPMS variable offset generates the same volume with CadQuery
    and VTK."""

    def variable_offset(x: np.ndarray, _: np.ndarray, __: np.ndarray):
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
    """Test if the TPMS variable offset raises an error when it is out of bounds."""

    def variable_offset(x: np.ndarray, _: np.ndarray, __: np.ndarray):
        return x + param  # x ∈ [-0.5, 0.5]

    # offset must be in [0, 2 * max(gyroid)]
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=variable_offset,
    )

    with pytest.raises((ValueError, NotImplementedError)):
        tpms.generate(type_part="sheet")


def test_tpms_generate_given_wrong_type_part_parameter_must_raise_ValueError():
    """Test if the TPMS type part raises an error when it is not recognized."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=0.5,
    )
    with pytest.raises(ValueError):
        tpms.generateVtk(type_part="fake")
    with pytest.raises(ValueError):
        tpms.generate(type_part="fake")


def test_tpms_given_wrong_cell_size_parameter_must_raise_ValueError():
    """Test if the TPMS cell size raises an error when it is out of bounds."""
    with pytest.raises(ValueError):
        microgen.Tpms(
            surface_function=microgen.surface_functions.gyroid,
            cell_size=(1, 1),
        )


def test_tpms_given_wrong_repeat_cell_parameter_must_raise_ValueError():
    """Test if the TPMS repeat cell raises an error when it is out of bounds."""
    with pytest.raises(ValueError):
        microgen.Tpms(
            surface_function=microgen.surface_functions.gyroid,
            repeat_cell=(1, 1, 1, 1),
        )


def test_tpms_given_wrong_density_parameter_must_raise_ValueError():
    """Test if the TPMS density raises an error when it is out of bounds."""
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
    """Test if the TPMS density generates the corresponding volume."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        density=0.2,
    )

    part = tpms.generateVtk(type_part=type_part)
    assert np.isclose(part.volume, tpms.grid.volume * 0.2, rtol=1e-2)


@pytest.mark.parametrize("part_type", ["lower skeletal", "upper skeletal", "sheet"])
def test_tpms_given_max_density_must_return_corresponding_offset(
    part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
):
    """Test if the TPMS max density returns the corresponding offset."""
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
    """Test if the TPMS with 100% sheet density returns a cube."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        density=1.0,
    )

    assert np.isclose(tpms.sheet.volume, tpms.grid.volume, rtol=1.0e-9)


def test_tpms_given_density_must_return_corresponding_offset():
    """Test if the TPMS density returns the corresponding offset."""
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
    """Test if the TPMS properties return the same value."""
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
    """Test if the TPMS surface is not empty."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
    )

    assert np.any(tpms.surface.points)
    assert np.any(tpms.surface.faces)


def test_tpms_given_negative_offset_for_skeletal_must_work_with_vtk_and_raise_error_with_cadquery():
    """Test if the TPMS with a negative offset generates a skeletal with VTK
    and raises an error with CadQuery."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=-1.0,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate(type_part="lower skeletal")

    sheet = tpms.generateVtk(type_part="lower skeletal").extract_surface()
    assert 0.0 < sheet.volume < np.abs(tpms.grid.volume)

    def including_negative_values(x: np.ndarray, _: np.ndarray, __: np.ndarray):
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
    """Test if the TPMS with a negative offset generates a sheet with VTK and
    raises an error with CadQuery."""

    def all_negative(x: np.ndarray, _: np.ndarray, __: np.ndarray):
        return -1.0 + x

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=all_negative,
    )
    with pytest.raises(ValueError):
        tpms.generate(type_part="sheet")

    assert tpms.generateVtk(type_part="sheet").volume == 0.0

    def including_negative_values(x: np.ndarray, _: np.ndarray, __: np.ndarray):
        return x

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=including_negative_values,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate(type_part="sheet")

    sheet = tpms.generateVtk(type_part="sheet").extract_surface()
    assert 0.0 < sheet.volume < np.abs(tpms.grid.volume)


def test_tpms_center_and_orientation_must_correspond():
    """Test if the center and orientation of the TPMS are correctly set."""
    center = (1.0, -2.0, 3.0)
    orientation = (np.pi / 3, -np.pi / 4, np.pi / 5)

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=0.5,
        center=center,
        orientation=orientation,
    )
    vtk_sheet = tpms.generateVtk(type_part="sheet")
    cad_sheet = tpms.generate(type_part="sheet")

    no_orientation = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=0.5,
        center=center,
    )

    assert np.allclose(vtk_sheet.center, center)
    assert np.allclose(cad_sheet.Center().toTuple(), center, rtol=1e-3)

    bbox = cad_sheet.BoundingBox()
    bounds = (bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax, bbox.zmin, bbox.zmax)
    assert np.allclose(bounds, vtk_sheet.bounds)

    assert not np.allclose(
        vtk_sheet.bounds,
        no_orientation.generateVtk(type_part="sheet").bounds,
    )
