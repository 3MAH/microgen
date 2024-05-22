"""Tests for the TPMS shapes generation."""

from __future__ import annotations

from inspect import getmembers, isfunction
from typing import Literal

import microgen
import numpy as np
import numpy.typing as npt
import pytest

TEST_DEFAULT_OFFSET = 0.5


def _get_microgen_surface_functions() -> list[str]:
    # Dont take into account deprecated surface functions named in camelCase
    return [
        func[0]
        for func in getmembers(microgen.surface_functions, isfunction)
        if not any(ele.isupper() for ele in func[0])
    ]


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal", "sheet"])
def test_tpms_given_cadquery_vtk_shapes_volume_must_be_equivalent(
    type_part: Literal["sheet", "lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    # Arrange
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        density=0.3,
    )

    # Act
    shape_cadquery = tpms.generate(type_part=type_part)
    shape_vtk = tpms.generate_vtk(type_part=type_part)

    # Assert
    assert np.isclose(shape_cadquery.Volume(), np.abs(shape_vtk.volume), rtol=1e-2)


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal"])
def test_tpms_given_cadquery_vtk_zero_offset_skeletals_volume_must_be_equivalent(
    type_part: Literal["lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS skeletals generated with CadQuery and VTK."""
    # Arrange
    tpms = microgen.Tpms(surface_function=microgen.surface_functions.schwarzP, offset=0)

    # Act
    shape_cadquery = tpms.generate(type_part=type_part)
    shape_vtk = tpms.generate_vtk(type_part=type_part)

    # Assert
    assert np.isclose(shape_cadquery.Volume(), np.abs(shape_vtk.volume), rtol=1e-2)


def test_tpms_given_non_default_cell_size_and_repeat_cell_must_have_same_volume_with_cad_and_vtk() -> (
    None
):
    """Test for non-default cell size and repeat cell values."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
        cell_size=(0.5, 2.0, 1.25),
        repeat_cell=(2, 1, 2),
    )

    shape_cadquery = tpms.generate(type_part="sheet")
    shape_vtk = tpms.generate_vtk(type_part="sheet")

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
    sheet = tpms.generate_vtk(type_part="sheet")
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
        surface_function=microgen.surface_functions.gyroid,
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
        surface_function=microgen.surface_functions.gyroid,
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
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
        repeat_cell=repeat_cell_zero,
    )

    tpms_repeat_max = coord_sys_tpms(
        radius=1.0,
        surface_function=microgen.surface_functions.gyroid,
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
    tpms = microgen.Tpms(surface_function=microgen.surface_functions.gyroid, offset=0)

    surface = tpms.generate(type_part="surface")
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
        surface_function=microgen.surface_functions.gyroid,
        offset=variable_offset,
    )

    shape_cadquery = tpms.generate(type_part="sheet", smoothing=0, verbose=True)
    shape_vtk = tpms.generate_vtk(type_part="sheet")

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
        surface_function=microgen.surface_functions.gyroid,
        offset=variable_offset,
    )

    with pytest.raises((ValueError, NotImplementedError)):
        tpms.generate(type_part="sheet")


def test_tpms_generate_given_wrong_type_part_parameter_must_raise_error() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
    )
    with pytest.raises(ValueError):
        tpms.generate_vtk(type_part="fake")
    with pytest.raises(ValueError):
        tpms.generate(type_part="fake")


def test_tpms_given_wrong_cell_size_parameter_must_raise_error() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    with pytest.raises(ValueError):
        microgen.Tpms(
            surface_function=microgen.surface_functions.gyroid,
            cell_size=(1, 1),
        )


def test_tpms_given_wrong_repeat_cell_parameter_must_raise_error() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    with pytest.raises(ValueError):
        microgen.Tpms(
            surface_function=microgen.surface_functions.gyroid,
            repeat_cell=(1, 1, 1, 1),
        )


def test_tpms_given_wrong_density_parameter_must_raise_error() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    with pytest.raises(ValueError):
        microgen.Tpms(
            surface_function=microgen.surface_functions.gyroid,
            density=0.0,
        )


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal", "sheet"])
def test_tpms_given_density_must_generate_tpms_with_correct_volume(
    type_part: Literal["sheet", "lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    expected_density = 0.2
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        density=expected_density,
    )

    part = tpms.generate_vtk(type_part=type_part)
    assert np.isclose(part.volume, abs(tpms.grid.volume) * expected_density, rtol=1e-2)


@pytest.mark.parametrize("type_part", ["lower skeletal", "upper skeletal", "sheet"])
def test_tpms_given_100_percent_density_must_return_a_cube(
    type_part: Literal["sheet", "lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        density=1.0,
    )

    assert np.isclose(
        tpms.generate_vtk(type_part=type_part).volume,
        abs(tpms.grid.volume),
        rtol=1.0e-9,
    )


def test_tpms_offset_from_density_given_density_must_return_valid_offset() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(surface_function=microgen.surface_functions.gyroid, offset=0)

    offset = microgen.Tpms.offset_from_density(
        surface_function=microgen.surface_functions.gyroid,
        density=0.5,
        part_type="sheet",
    )

    max_offset = -2.0 * np.min(tpms.grid["surface"])
    assert 0 < offset < max_offset


def test_tpms_given_property_must_return_the_same_value() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        density=0.2,
    )
    skeletals = tpms.skeletals

    assert tpms.upper_skeletal == skeletals[0]
    assert tpms.lower_skeletal == skeletals[1]
    assert tpms.sheet == tpms.sheet
    assert tpms.generate_vtk(type_part="surface") == tpms.surface


def test_tpms_given_surface_must_not_be_empty() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(surface_function=microgen.surface_functions.gyroid, offset=0)

    assert np.any(tpms.surface.points)
    assert np.any(tpms.surface.faces)


def test_tpms_given_negative_offset_for_skeletal_must_work_with_vtk_and_raise_error_with_cadquery() -> (
    None
):
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=-1.0,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate(type_part="lower skeletal")

    sheet = tpms.generate_vtk(type_part="lower skeletal")
    assert 0.0 < sheet.volume < np.abs(tpms.grid.volume)

    def including_negative_values(
        x: npt.NDArray[np.float64],
        _: npt.NDArray[np.float64],
        __: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return x

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=including_negative_values,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate(type_part="lower skeletal")

    sheet = tpms.generate_vtk(type_part="lower skeletal")
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
        surface_function=microgen.surface_functions.gyroid,
        offset=all_negative,
    )
    with pytest.raises(ValueError):
        tpms.generate(type_part="sheet")

    assert tpms.generate_vtk(type_part="sheet").volume == 0.0

    def including_negative_values(
        x: npt.NDArray[np.float64],
        _: npt.NDArray[np.float64],
        __: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return x

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=including_negative_values,
    )
    with pytest.raises(NotImplementedError):
        tpms.generate(type_part="sheet")

    sheet = tpms.generate_vtk(type_part="sheet")
    assert 0.0 < sheet.volume < np.abs(tpms.grid.volume)


def test_tpms_center_and_orientation_must_correspond() -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    center = (1.0, -2.0, 3.0)
    orientation = (np.pi / 3, -np.pi / 4, np.pi / 5)

    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
        center=center,
        orientation=orientation,
    )
    vtk_sheet = tpms.generate_vtk(type_part="sheet")
    cad_sheet = tpms.generate(type_part="sheet")

    no_orientation = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
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
        no_orientation.generate_vtk(type_part="sheet").bounds,
    )


@pytest.mark.parametrize("part_type", ["sheet", "lower skeletal", "upper skeletal"])
def test_tpms_generate_vtk_check_that_volume_has_changed_when_the_offset_is_updated(
    part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
    )
    first_part = tpms.generate_vtk(type_part=part_type)

    tpms.offset *= 2.0
    second_part = tpms.generate_vtk(type_part=part_type)
    assert not np.isclose(first_part.volume, second_part.volume)


@pytest.mark.parametrize("part_type", ["sheet", "lower skeletal", "upper skeletal"])
def test_tpms_generate_grid_vtk_check_that_volume_has_changed_when_the_offset_is_updated(
    part_type: Literal["sheet", "lower skeletal", "upper skeletal"],
) -> None:
    """Test for the volume of the TPMS shapes generated with CadQuery and VTK."""
    tpms = microgen.Tpms(
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
    )
    first_part = tpms.generate_grid_vtk(type_part=part_type)

    tpms.offset *= 2.0
    second_part = tpms.generate_grid_vtk(type_part=part_type)
    assert not np.isclose(first_part.volume, second_part.volume)


def test_infill_given_cell_size_must_use_corresponding_repeat_cell() -> None:
    """Test if the repeat cell is computed correctly."""
    tpms = microgen.Infill(
        obj=microgen.Box(dim=(1.0, 1.0, 1.0)).generate_vtk(),
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
        cell_size=(0.5, 1.0, 1.0),
    )
    expected_repeat_cell = (2, 1, 1)
    assert np.allclose(tpms.repeat_cell, expected_repeat_cell)


def test_infill_given_repeat_cell_must_use_corresponding_cell_size() -> None:
    """Test if the cell size is computed correctly."""
    tpms = microgen.Infill(
        obj=microgen.Box(dim=(1.0, 1.0, 1.0)).generate_vtk(),
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
        repeat_cell=(1, 1, 2),
    )

    expected_cell_size = (1.0, 1.0, 0.5)
    assert np.allclose(tpms.cell_size, expected_cell_size)


@pytest.mark.parametrize("kwarg", [{"cell_size": 0.5}, {"repeat_cell": 2}])
def test_infill_bounds_match_obj_bounds(kwarg: dict[str, int | float]) -> None:
    """Test if the grid bounds match the object bounds."""
    obj = microgen.Ellipsoid(radii=(1.0, 2.0 / 3.0, 0.5)).generate_vtk()
    tpms = microgen.Infill(
        obj=obj,
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
        **kwarg,
    )

    grid_bounds = np.array(tpms.grid.bounds)
    grid_dim = grid_bounds[1::2] - grid_bounds[::2]

    obj_bounds = np.array(obj.bounds)
    obj_dim = obj_bounds[1::2] - obj_bounds[::2]

    assert np.all(obj_dim <= grid_dim) or np.allclose(obj_dim, grid_dim, rtol=1e-2)
    assert np.all(grid_dim < obj_dim + tpms.cell_size)


def test_infill_fills_the_object_with_any_normal_orientation() -> None:
    """Test if the object is filled correctly with any normal orientation."""
    mesh = microgen.Box(dim=(1.0, 1.0, 1.0)).generate_vtk()
    first = microgen.Infill(
        obj=mesh,
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
        cell_size=0.5,
    )
    mesh = mesh.triangulate()
    mesh.flip_normals()
    second = microgen.Infill(
        obj=mesh,
        surface_function=microgen.surface_functions.gyroid,
        offset=TEST_DEFAULT_OFFSET,
        cell_size=0.5,
    )

    assert np.allclose(first.sheet.volume, second.sheet.volume)
    assert first.sheet.points.shape == second.sheet.points.shape


def test_infill_given_repeat_cell_and_cell_size_must_raise_an_error() -> None:
    """Test if the cell size is computed correctly."""
    with pytest.raises(ValueError):
        microgen.Infill(
            obj=microgen.Box(dim=(1.0, 1.0, 1.0)).generate_vtk(),
            surface_function=microgen.surface_functions.gyroid,
            offset=TEST_DEFAULT_OFFSET,
            repeat_cell=(2, 1, 2),
            cell_size=(0.5, 1.0, 0.5),
        )


def test_infill_raises_error_when_cell_size_is_too_large() -> None:
    """Test if the cell size is too large compared to the given object size."""
    too_large_cell_size = (1.0, 2.0, 1.0)
    with pytest.raises(ValueError):
        microgen.Infill(
            obj=microgen.Box(dim=(1.0, 1.0, 1.0)).generate_vtk(),
            surface_function=microgen.surface_functions.gyroid,
            offset=TEST_DEFAULT_OFFSET,
            cell_size=too_large_cell_size,
        )


def test_tpms_both_offset_and_density_given_must_raise_error() -> None:
    """Test whether providing both offset and density results in an error."""
    with pytest.raises(ValueError):
        microgen.Tpms(
            surface_function=microgen.surface_functions.gyroid,
            offset=TEST_DEFAULT_OFFSET,
            density=0.5,
        )


def test_tpms_none_offset_and_density_given_must_raise_error() -> None:
    """Test whether omitting both offset and density results in an error."""
    with pytest.raises(ValueError):
        microgen.Tpms(
            surface_function=microgen.surface_functions.gyroid,
        )
