"""Test the lattice shapes."""

from typing import Literal

import cadquery as cq
import numpy as np
import pytest

from microgen.shape.strut_lattice import (
    BodyCenteredCubic,
    Cubic,
    Cuboctahedron,
    Diamond,
    FaceCenteredCubic,
    Octahedron,
    OctetTruss,
    RhombicCuboctahedron,
    RhombicDodecahedron,
    TruncatedCube,
    TruncatedCuboctahedron,
    TruncatedOctahedron,
)

from .test_lattice_utils import PRESET_LATTICE_DATA


@pytest.mark.parametrize(
    "shape",
    [
        BodyCenteredCubic,
        Cubic,
        Cuboctahedron,
        Diamond,
        FaceCenteredCubic,
        Octahedron,
        OctetTruss,
        RhombicCuboctahedron,
        RhombicDodecahedron,
        TruncatedCube,
        TruncatedCuboctahedron,
        TruncatedOctahedron,
    ],
)
def test_lattice_vertices_strut_centers_and_directions_must_correspond_to_preset_lattice_data(
    shape,
) -> None:
    lattice = shape(strut_radius=0.05, cell_size=1.0, strut_joints=False)

    assert np.allclose(
        np.sort(lattice.vertices.flat),
        np.sort(PRESET_LATTICE_DATA[shape.__name__].vertices.flat),
    )
    assert np.allclose(
        np.sort(lattice.strut_centers.flat),
        np.sort(PRESET_LATTICE_DATA[shape.__name__].strut_centers.flat),
    )
    assert np.allclose(
        np.sort(lattice.strut_directions_cartesian.flat),
        np.sort(PRESET_LATTICE_DATA[shape.__name__].strut_directions.flat),
    )


@pytest.mark.parametrize("shape", [Diamond, Cubic, OctetTruss])
@pytest.mark.parametrize("cell_size", [0.5, 1.0, 2.0])
@pytest.mark.parametrize("expected_density", [0.25, 0.8])
def test_lattice_given_density_and_cell_size_must_match_computed_density(
    shape: Literal[Diamond, Cubic, OctetTruss],
    expected_density: float,
    cell_size: float,
) -> None:
    """Test for the density of lattice shapes generated with CadQuery and VTK."""
    # Arrange
    shape_fit_to_density = shape(density=expected_density, cell_size=cell_size)

    # Act
    shape_density = shape_fit_to_density.density

    # Assert
    assert np.isclose(expected_density, shape_density, rtol=0.01)
