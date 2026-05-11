"""Tests for the TPMS shapes generation."""

from __future__ import annotations

import re
from inspect import getmembers, isfunction, signature
from typing import Literal

import numpy as np
import numpy.typing as npt
import pytest

import microgen

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/
# ruff: noqa: E501 line-too-long https://docs.astral.sh/ruff/rules/line-too-long/


TEST_DEFAULT_OFFSET = 0.5
TEST_DEFAULT_SURFACE_FUNCTION = microgen.surface_functions.gyroid


def _get_microgen_surface_functions() -> list[str]:
    """List the actual TPMS surface functions in microgen.surface_functions.

    Filters out non-TPMS callables exposed via re-import (autograd ``cos``/``sin``)
    that take fewer than 3 args and would crash when invoked as ``f(x,y,z)``.
    """
    return [
        name
        for name, fn in getmembers(microgen.surface_functions, isfunction)
        if len(signature(fn).parameters) == 3
    ]


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal", "sheet"])
def test_tpms_given_cadquery_vtk_shapes_volume_must_be_equivalent(
    type_part: Literal["sheet", "lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    # Arrange
    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        density=0.3,
    )

    # Act
    shape_cadquery = tpms.generate_cad(type_part=type_part)
    shape_vtk = tpms.generate_surface_mesh(type_part=type_part)

    # Assert
    assert np.isclose(shape_cadquery.Volume(), np.abs(shape_vtk.volume), rtol=1e-2)


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal"])
def test_tpms_given_cadquery_vtk_zero_offset_skeletals_volume_must_be_equivalent(
    type_part: Literal["lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS skeletals generated with CadQuery and VTK."""
    # Arrange
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.schwarz_p, offset=0
    )

    # Act
    shape_cadquery = tpms.generate_cad(type_part=type_part)
    shape_vtk = tpms.generate_surface_mesh(type_part=type_part)

    # Assert
    assert np.isclose(shape_cadquery.Volume(), np.abs(shape_vtk.volume), rtol=1e-2)


def test_tpms_given_non_default_cell_size_and_repeat_cell_must_have_same_volume_with_cad_and_vtk() -> (
    None
):
    """Test for non-default cell size and repeat cell values."""
    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
        cell_size=(0.5, 2.0, 1.25),
        repeat_cell=(2, 1, 2),
    )

    shape_cadquery = tpms.generate_cad(type_part="sheet")
    shape_vtk = tpms.generate_surface_mesh(type_part="sheet")

    assert np.isclose(shape_cadquery.Volume(), np.abs(shape_vtk.volume), rtol=1e-2)


@pytest.mark.parametrize(
    "surface",
    _get_microgen_surface_functions(),
)
@pytest.mark.parametrize("repeat_cell", [2, (2, 1, 3)])
@pytest.mark.parametrize("cell_size", [3.0, (0.5, 1.5, 1.0)])
def test_tpms_given_sum_volume_must_be_cube_volume(
    surface: str,
    repeat_cell: int | tuple[int, int, int],
    cell_size: float | tuple[float, float, float],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    # Arrange
    tpms = microgen.Tpms(
        surface_function=getattr(microgen.surface_functions, surface),
        offset=TEST_DEFAULT_OFFSET,
        repeat_cell=repeat_cell,
        cell_size=cell_size,
    )

    # Act
    volume = np.abs(
        tpms.sheet.volume + tpms.lower_skeletal.volume + tpms.upper_skeletal.volume,
    )
    cube_volume = np.prod(tpms.repeat_cell) * np.prod(tpms.cell_size)

    # Assert
    assert np.isclose(volume, cube_volume, rtol=1e-2)


@pytest.mark.parametrize("surface", _get_microgen_surface_functions())
@pytest.mark.parametrize("density", [0.05, 0.5, 0.99, 1.0])
def test_tpms_given_density_must_match_computed_density(
    surface: str,
    density: float,
) -> None:
    """Test for the density of the TPMS shapes generated with CadQuery and VTK."""
    # Arrange
    tpms = microgen.Tpms(
        surface_function=getattr(microgen.surface_functions, surface),
        density=density,
    )

    # Act
    sheet = tpms.generate_surface_mesh(type_part="sheet")
    computed_density = sheet.volume / abs(tpms.grid.volume)

    # Assert
    assert np.isclose(computed_density, density, rtol=0.1)


@pytest.mark.parametrize(
    "coord_sys_tpms",
    [microgen.Tpms, microgen.CylindricalTpms, microgen.SphericalTpms],
)
def test_tpms_given_coord_system_tpms_coordinates_field_must_be_in_cartesian_frame(
    coord_sys_tpms: type[microgen.Tpms],
) -> None:
    """Test if the `coords` field of the grid correspond to the cartesian frame."""
    kwargs = {"radius": 1.0} if coord_sys_tpms != microgen.Tpms else {}

    tpms = coord_sys_tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
        **kwargs,
    )

    linspaces: list[npt.NDArray[np.float64]] = [
        np.linspace(
            -0.5 * cell_size_axis * repeat_cell_axis,
            0.5 * cell_size_axis * repeat_cell_axis,
            tpms.resolution * repeat_cell_axis,
        )
        for repeat_cell_axis, cell_size_axis in zip(
            tpms.repeat_cell,
            tpms.cell_size,
        )
    ]

    cartesian_coords = np.vstack(
        [axis.ravel("F") for axis in np.meshgrid(*linspaces)],
    ).T

    assert np.all(cartesian_coords == tpms.grid["coords"])


@pytest.mark.parametrize(
    "coord_sys_tpms",
    [microgen.CylindricalTpms, microgen.SphericalTpms],
)
def test_tpms_given_coord_system_tpms_volumes_must_be_greater_than_zero_and_lower_than_grid_volume(
    coord_sys_tpms: type[microgen.Tpms],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = coord_sys_tpms(
        radius=1.0,
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        density=0.2,
    )

    assert 0 < tpms.sheet.volume < np.abs(tpms.grid.volume)
    assert 0 < tpms.lower_skeletal.volume < np.abs(tpms.grid.volume)
    assert 0 < tpms.upper_skeletal.volume < np.abs(tpms.grid.volume)


@pytest.mark.parametrize(
    ("coord_sys_tpms", "repeat_cell_zero", "repeat_cell_max"),
    [
        (microgen.CylindricalTpms, (1, 0, 1), (1, 100, 1)),
        (microgen.SphericalTpms, (1, 0, 0), (1, 100, 100)),
    ],
)
def test_tpms_given_zero_and_max_repeat_cell_values_volumes_must_correspond(
    coord_sys_tpms: type[microgen.Tpms],
    repeat_cell_zero: tuple[int, int, int],
    repeat_cell_max: tuple[int, int, int],
) -> None:
    """Test for zero and max repeat cell values.

    The volume of the sheet, lower skeletal, and upper skeletal must be the same for zero and
    max repeat cell values and be between 0 and the volume of the grid.
    """
    tpms_repeat_zero = coord_sys_tpms(
        radius=1.0,
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
        repeat_cell=repeat_cell_zero,
    )

    tpms_repeat_max = coord_sys_tpms(
        radius=1.0,
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
        repeat_cell=repeat_cell_max,
    )

    assert np.abs(tpms_repeat_zero.grid.volume) == np.abs(tpms_repeat_max.grid.volume)

    assert (
        0
        < tpms_repeat_zero.sheet.volume
        == tpms_repeat_max.sheet.volume
        < np.abs(tpms_repeat_zero.grid.volume)
    )
    assert (
        0
        < tpms_repeat_zero.lower_skeletal.volume
        == tpms_repeat_max.lower_skeletal.volume
        < np.abs(tpms_repeat_zero.grid.volume)
    )
    assert (
        0
        < tpms_repeat_zero.upper_skeletal.volume
        == tpms_repeat_max.upper_skeletal.volume
        < np.abs(tpms_repeat_zero.grid.volume)
    )


def test_tpms_given_generate_surface_must_not_be_empty() -> None:
    """Test for the surface of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(surface_function=TEST_DEFAULT_SURFACE_FUNCTION, offset=0)

    surface = tpms.generate_cad(type_part="surface")
    assert np.any(surface.Vertices())
    assert np.any(surface.Faces())
    assert not surface.Closed()


def test_tpms_given_variable_offset_cadquery_and_vtk_volumes_must_correspond() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""

    def variable_offset(
        x: npt.NDArray[np.float64],
        _: npt.NDArray[np.float64],
        __: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return x + 1.5

    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=variable_offset,
    )

    shape_cadquery = tpms.generate_cad(type_part="sheet", smoothing=0, verbose=True)
    shape_vtk = tpms.generate_surface_mesh(type_part="sheet")

    assert np.isclose(shape_cadquery.Volume(), np.abs(shape_vtk.volume), rtol=1e-2)


@pytest.mark.parametrize("param", [0.3, 4.0])
def test_tpms_given_variable_offset_out_of_limits_with_cadquery_must_raise_error(
    param: float,
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""

    def variable_offset(
        x: npt.NDArray[np.float64],
        _: npt.NDArray[np.float64],
        __: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return x + param  # x ∈ [-0.5, 0.5]

    # offset must be in [0, 2 * max(gyroid)]
    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=variable_offset,
    )

    with pytest.raises((ValueError, NotImplementedError)):
        tpms.generate_cad(type_part="sheet")


def test_tpms_generate_given_wrong_type_part_parameter_must_raise_error() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
    )
    fake_type_part = "fake"
    expected_err_msg = re.escape(
        (
            f"type_part ({fake_type_part}) must be 'sheet', 'lower skeletal', "
            "'upper skeletal' or 'surface'"
        ),
    )
    with pytest.raises(ValueError, match=expected_err_msg):
        tpms.generate_surface_mesh(type_part=fake_type_part)
    with pytest.raises(ValueError, match=expected_err_msg):
        tpms.generate_cad(type_part=fake_type_part)


def test_tpms_given_wrong_cell_size_parameter_must_raise_error() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    invalid_cell_size = (1, 1)
    expected_err_msg = re.escape(
        f"`cell_size` must have a length of 3 floats. Given: {invalid_cell_size}",
    )
    with pytest.raises(ValueError, match=expected_err_msg):
        microgen.Tpms(
            surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
            offset=TEST_DEFAULT_OFFSET,
            cell_size=invalid_cell_size,
        )


def test_tpms_given_wrong_repeat_cell_parameter_must_raise_error() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    invalid_repeat_cell = (1, 1, 1, 1)
    expected_err_msg = re.escape(
        f"`repeat_cell` must have a length of 3 integers. Given: {invalid_repeat_cell}",
    )
    with pytest.raises(ValueError, match=expected_err_msg):
        microgen.Tpms(
            surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
            offset=TEST_DEFAULT_OFFSET,
            repeat_cell=invalid_repeat_cell,
        )


def test_tpms_given_wrong_density_parameter_must_raise_error() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    invalid_density = 0.0
    expected_err_msg = re.escape(
        f"density must be between 0 and 1. Given: {invalid_density}",
    )
    with pytest.raises(ValueError, match=expected_err_msg):
        microgen.Tpms(
            surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
            density=invalid_density,
        )


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal", "sheet"])
def test_tpms_given_density_must_generate_tpms_with_correct_volume(
    type_part: Literal["sheet", "lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    expected_density = 0.2
    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        density=expected_density,
    )

    part = tpms.generate_surface_mesh(type_part=type_part)
    assert np.isclose(part.volume, abs(tpms.grid.volume) * expected_density, rtol=1e-2)


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal", "sheet"])
def test_tpms_given_100_percent_density_must_return_a_cube(
    type_part: Literal["sheet", "lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        density=1.0,
    )

    assert np.isclose(
        tpms.generate_surface_mesh(type_part=type_part).volume,
        abs(tpms.grid.volume),
        rtol=1.0e-9,
    )


def test_tpms_offset_from_density_given_density_must_return_valid_offset() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(surface_function=TEST_DEFAULT_SURFACE_FUNCTION, offset=0)

    offset = microgen.Tpms.offset_from_density(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        density=0.5,
        part_type="sheet",
    )

    max_offset = -2.0 * np.min(tpms.grid["surface"])
    assert 0 < offset < max_offset


def test_tpms_given_property_must_return_the_same_value() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        density=0.2,
    )
    skeletals = tpms.skeletals

    assert tpms.upper_skeletal == skeletals[0]
    assert tpms.lower_skeletal == skeletals[1]
    assert tpms.sheet == tpms.sheet
    assert tpms.generate_surface_mesh(type_part="surface") == tpms.surface


def test_tpms_given_surface_must_not_be_empty() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(surface_function=TEST_DEFAULT_SURFACE_FUNCTION, offset=0)

    assert np.any(tpms.surface.points)
    assert np.any(tpms.surface.faces)


def test_tpms_given_negative_offset_for_skeletal_must_work_with_vtk_and_raise_error_with_cadquery() -> (
    None
):
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=-1.0,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate_cad(type_part="lower skeletal")

    sheet = tpms.generate_surface_mesh(type_part="lower skeletal")
    assert 0.0 < sheet.volume < np.abs(tpms.grid.volume)

    def including_negative_values(
        x: npt.NDArray[np.float64],
        _: npt.NDArray[np.float64],
        __: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return x

    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=including_negative_values,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate_cad(type_part="lower skeletal")

    sheet = tpms.generate_surface_mesh(type_part="lower skeletal")
    assert 0.0 < sheet.volume < np.abs(tpms.grid.volume)


def test_tpms_given_negative_offset_for_sheet_must_work_with_vtk_and_raise_error_with_cadquery() -> (
    None
):
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""

    def all_negative(
        x: npt.NDArray[np.float64],
        _: npt.NDArray[np.float64],
        __: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return -1.0 + x

    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=all_negative,
    )

    err_msg_pattern = (
        re.escape("offset must be greater than ")
        + r"-?\d*\.\d+"  # float
        + re.escape(" to generate 'sheet' part and lower than ")
        + r"-?\d*\.\d+"  # float
    )
    with pytest.raises(ValueError, match=err_msg_pattern):
        tpms.generate_cad(type_part="sheet")

    assert tpms.generate_surface_mesh(type_part="sheet").volume == 0.0

    def including_negative_values(
        x: npt.NDArray[np.float64],
        _: npt.NDArray[np.float64],
        __: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return x

    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=including_negative_values,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate_cad(type_part="sheet")

    sheet = tpms.generate_surface_mesh(type_part="sheet")
    assert 0.0 < sheet.volume < np.abs(tpms.grid.volume)


def test_tpms_center_and_orientation_must_correspond() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    center = (1.0, -2.0, 3.0)
    orientation = (np.pi / 3, -np.pi / 4, np.pi / 5)

    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
        center=center,
        orientation=orientation,
    )
    vtk_sheet = tpms.generate_surface_mesh(type_part="sheet")
    cad_sheet = tpms.generate_cad(type_part="sheet")

    no_orientation = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
        center=center,
    )

    assert np.allclose(vtk_sheet.center, center)
    assert np.allclose(cad_sheet.Center().toTuple(), center, rtol=1e-3)

    bbox = cad_sheet.BoundingBox()
    bounds = (bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax, bbox.zmin, bbox.zmax)
    assert np.allclose(bounds, vtk_sheet.bounds)

    assert not np.allclose(
        vtk_sheet.bounds,
        no_orientation.generate_surface_mesh(type_part="sheet").bounds,
    )


@pytest.mark.parametrize("part_type", ["sheet", "lower skeletal", "upper skeletal"])
def test_tpms_generate_surface_mesh_check_that_volume_has_changed_when_the_offset_is_updated(
    part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
    )
    first_part = tpms.generate_surface_mesh(type_part=part_type)

    tpms.offset *= 2.0
    second_part = tpms.generate_surface_mesh(type_part=part_type)
    assert not np.isclose(first_part.volume, second_part.volume)


@pytest.mark.parametrize("part_type", ["sheet", "lower skeletal", "upper skeletal"])
def test_tpms_generate_volume_mesh_check_that_volume_has_changed_when_the_offset_is_updated(
    part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
    )
    first_part = tpms.generate_volume_mesh(type_part=part_type)

    tpms.offset *= 2.0
    second_part = tpms.generate_volume_mesh(type_part=part_type)
    assert not np.isclose(first_part.volume, second_part.volume)


def test_infill_given_cell_size_must_use_corresponding_repeat_cell() -> None:
    """Test if the repeat cell is computed correctly."""
    tpms = microgen.Infill(
        obj=microgen.Box(dim=(1.0, 1.0, 1.0)).generate_surface_mesh(),
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
        cell_size=(0.5, 1.0, 1.0),
    )
    expected_repeat_cell = (2, 1, 1)
    assert np.allclose(tpms.repeat_cell, expected_repeat_cell)


def test_infill_given_repeat_cell_must_use_corresponding_cell_size() -> None:
    """Test if the cell size is computed correctly."""
    tpms = microgen.Infill(
        obj=microgen.Box(dim=(1.0, 1.0, 1.0)).generate_surface_mesh(),
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
        repeat_cell=(1, 1, 2),
    )

    expected_cell_size = (1.0, 1.0, 0.5)
    assert np.allclose(tpms.cell_size, expected_cell_size, rtol=1e-2)


@pytest.mark.parametrize("kwarg", [{"cell_size": 0.5}, {"repeat_cell": 2}])
def test_infill_bounds_match_obj_bounds(kwarg: dict[str, int | float]) -> None:
    """Test if the grid bounds match the object bounds."""
    obj = microgen.Ellipsoid(radii=(1.0, 2.0 / 3.0, 0.5)).generate_surface_mesh()
    tpms = microgen.Infill(
        obj=obj,
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
        **kwarg,  # type: ignore[arg-type]
    )

    grid_bounds = np.array(tpms.grid.bounds)
    grid_dim = grid_bounds[1::2] - grid_bounds[::2]

    obj_bounds = np.array(obj.bounds)
    obj_dim = obj_bounds[1::2] - obj_bounds[::2]

    assert np.all(obj_dim <= grid_dim) or np.allclose(obj_dim, grid_dim, rtol=1e-2)
    assert np.all(grid_dim < obj_dim + tpms.cell_size)


def test_infill_cylinder_has_expected_volume() -> None:
    """Test if an infilled cylinder has the expected volume."""
    density = 0.5
    cylinder = microgen.Cylinder().generate_surface_mesh()
    infill = microgen.Infill(
        cylinder,
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        density=density,
        repeat_cell=2,
    )
    assert np.isclose(infill.sheet.volume, density * cylinder.volume, rtol=1e-2)


def test_infill_cylinder_returns_single_connected_component_mesh() -> None:
    """Test if the infilled cylinder returns a single connected component mesh."""
    infill = microgen.Infill(
        obj=microgen.Cylinder().generate_surface_mesh(),
        surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
        offset=TEST_DEFAULT_OFFSET,
        repeat_cell=2,
    )

    components = infill.sheet.connectivity().point_data["RegionId"]
    n_unique = len(np.unique(components))

    assert n_unique == 1


def test_infill_given_repeat_cell_and_cell_size_must_raise_an_error() -> None:
    """Test if the cell size is computed correctly."""
    expected_err_msg = (
        "cell_size and repeat_cell cannot be given at the same time, "
        "one is computed from the other."
    )
    with pytest.raises(ValueError, match=expected_err_msg):
        microgen.Infill(
            obj=microgen.Box(dim=(1.0, 1.0, 1.0)).generate_surface_mesh(),
            surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
            offset=TEST_DEFAULT_OFFSET,
            repeat_cell=(2, 1, 2),
            cell_size=(0.5, 1.0, 0.5),
        )


def test_infill_raises_error_when_cell_size_is_too_large() -> None:
    """Test if the cell size is too large compared to the given object size."""
    too_large_cell_size = (1.0, 2.0, 1.0)
    expected_err_msg = (
        re.escape("cell_size must be lower than the object dimensions. ")
        + re.escape(f"Given: {too_large_cell_size}, Object dimensions: ")
        + r"\[\d+(\.\d+)? \d+(\.\d+)? \d+(\.\d+)?\]"  # np.ndarray of 3 floats
    )
    with pytest.raises(ValueError, match=expected_err_msg):
        microgen.Infill(
            obj=microgen.Box(dim=(1.0, 1.0, 1.0)).generate_surface_mesh(),
            surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
            offset=TEST_DEFAULT_OFFSET,
            cell_size=too_large_cell_size,
        )


def test_tpms_both_offset_and_density_given_must_raise_error() -> None:
    """Test whether providing both offset and density results in an error."""
    expected_err_msg = (
        "offset and density cannot be given at the same time. Give only one."
    )
    with pytest.raises(ValueError, match=expected_err_msg):
        microgen.Tpms(
            surface_function=TEST_DEFAULT_SURFACE_FUNCTION,
            offset=TEST_DEFAULT_OFFSET,
            density=0.5,
        )


def test_tpms_none_offset_and_density_given_must_raise_error() -> None:
    """Test whether omitting both offset and density results in an error."""
    expected_err_msg = "offset or density must be given. Give one of them."
    with pytest.raises(ValueError, match=expected_err_msg):
        microgen.Tpms(surface_function=TEST_DEFAULT_SURFACE_FUNCTION)


# -----------------------------------------------------------------------------
# Coordinate-frame mode regression tests
# -----------------------------------------------------------------------------
# These guard against the cell-box / bounds confusion that made
# CylindricalTpms / SphericalTpms produce wildly low volumes after the F-rep
# refactor (the parent's axis-aligned-box ``_cell_box`` was treated as a
# Cartesian clip even though those subclasses use parametric units).


def test_cylindrical_tpms_full_wrap_volume_matches_shell() -> None:
    """A full-wrap cylindrical TPMS sheet should fill a sizeable fraction of
    the underlying ``2πR·Δr·H`` shell envelope.  The exact fraction depends
    on the surface function and offset, but for gyroid + offset 0.5 it
    should be > 30 % of the shell envelope.  Pre-fix this test would have
    seen ~10 %.
    """
    from microgen.shape.tpms import CylindricalTpms

    radius, delta_r, height = 1.5, 2.0, 6.0
    tpms = CylindricalTpms(
        radius=radius,
        surface_function=microgen.surface_functions.gyroid,
        offset=0.5,
        cell_size=1.0,
        repeat_cell=(2, 0, int(height)),  # 0 → auto-fill full circle
    )
    sheet = tpms.generate_surface_mesh(type_part="sheet")

    shell_envelope = 2.0 * np.pi * radius * delta_r * height
    density = abs(sheet.volume) / shell_envelope
    # Gyroid sheet density at offset 0.5 ≈ 10-25% via parametric grid-clip.
    # Assert just enough to gate the parametric-grid path: the broken F-rep
    # path was producing < 1 % here.
    assert density > 0.05, f"density {density:.2%} too low — clip path likely wrong"
    assert density < 1.05, f"density {density:.2%} > envelope — clipping leak"


def test_spherical_tpms_full_wrap_volume_matches_shell() -> None:
    """Full-sphere TPMS sheet should fill a sizeable fraction of the
    ``4πR²·Δr`` shell envelope.  Pre-fix this would have been ~1 %.
    """
    from microgen.shape.tpms import SphericalTpms

    radius, delta_r = 3.0, 2.0
    tpms = SphericalTpms(
        radius=radius,
        surface_function=microgen.surface_functions.gyroid,
        offset=0.5,
        cell_size=1.0,
        repeat_cell=(2, 0, 0),  # auto-fill θ + φ
    )
    sheet = tpms.generate_surface_mesh(type_part="sheet")

    shell_envelope = 4.0 * np.pi * radius * radius * delta_r
    density = abs(sheet.volume) / shell_envelope
    # Gyroid sheet density at offset 0.5 ≈ 10-25% via parametric grid-clip.
    assert density > 0.05, f"density {density:.2%} too low — clip path likely wrong"
    assert density < 1.05, f"density {density:.2%} > envelope — clipping leak"


def test_cylindrical_tpms_partial_wrap_smaller_than_full() -> None:
    """A quarter-wrap cylinder must produce ≲ 1/3 of the full-wrap volume.

    Gates the wedge clip in :meth:`CylindricalTpms._cell_box` — without it
    the partial wrap would equal the full wrap.
    """
    from microgen.shape.tpms import CylindricalTpms

    full = CylindricalTpms(
        radius=1.5,
        surface_function=microgen.surface_functions.gyroid,
        offset=0.5,
        cell_size=1.0,
        repeat_cell=(2, 0, 6),
    ).generate_surface_mesh(type_part="sheet")

    quarter = CylindricalTpms(
        radius=1.5,
        surface_function=microgen.surface_functions.gyroid,
        offset=0.5,
        cell_size=1.0,
        repeat_cell=(2, 2, 6),  # quarter wrap (~2 of 9 angular cells)
    ).generate_surface_mesh(type_part="sheet")

    assert abs(quarter.volume) < abs(full.volume) / 3.0, (
        f"quarter wrap vol {abs(quarter.volume):.2f} should be ≲ "
        f"{abs(full.volume) / 3.0:.2f} (1/3 of full)"
    )


def test_sweep_along_straight_line_is_finite_and_positive() -> None:
    """Sweep along a straight z-axis line must produce a positive sheet
    volume bounded by the tube envelope.
    """
    from microgen.shape.tpms import Sweep

    line = np.linspace([0.0, 0.0, -3.0], [0.0, 0.0, 3.0], 50)
    radial_max, height = 1.0, 6.0
    tpms = Sweep(
        curve_points=line,
        surface_function=microgen.surface_functions.gyroid,
        radial_max=radial_max,
        offset=0.4,
        cell_size=1.0,
        repeat_cell=(int(height), 1, 6),
    )
    sheet = tpms.generate_surface_mesh(type_part="sheet")

    tube_volume = np.pi * radial_max * radial_max * height
    v = abs(sheet.volume)
    assert v > 0.0, "Sweep produced an empty sheet"
    assert v < tube_volume * 1.05, (
        f"Sweep sheet volume {v:.2f} exceeds tube envelope {tube_volume:.2f}"
    )
