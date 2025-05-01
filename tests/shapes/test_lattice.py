"""Test the lattice shapes."""

import numpy as np
from typing import Literal
import pytest
import cadquery as cq

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
@pytest.mark.parametrize("density", [0.25, 0.8])
def test_lattice_diamond_given_density_and_cell_size_must_match_computed_density(
    shape: Literal[Diamond, Cubic, OctetTruss],
    density: float,
    cell_size: float,
) -> None:
    """Test for the density of the diamond lattice shape generated with CadQuery and VTK."""
    # Arrange
    expected_density = density

    # Act
    shape_fit_to_density = shape(density=expected_density, cell_size=cell_size)
    density = shape_fit_to_density.density

    generated_cad = shape_fit_to_density.cad_shape

    # Assert
    assert isinstance(generated_cad, cq.Shape)
    assert np.isclose(expected_density, density, rtol=0.01)
