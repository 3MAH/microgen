"""Test the lattice shapes."""

from typing import Literal

import cadquery as cq
import numpy as np
import pytest

from microgen import is_periodic
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


def test_lattice_generate_vtk_periodic_must_produce_periodic_mesh() -> None:
    """generate_vtk(periodic=True) must produce a mesh with periodic node positions."""
    lattice = OctetTruss(strut_radius=0.05, cell_size=1.0)
    mesh = lattice.generate_vtk(size=0.1, periodic=True)

    assert is_periodic(
        mesh.points
    ), "Mesh generated with periodic=True must be periodic"


def test_lattice_generate_vtk_must_not_reuse_non_periodic_mesh_for_periodic_request() -> None:
    """generate_vtk must cache by parameters, not by instance only."""
    lattice = OctetTruss(strut_radius=0.05, cell_size=1.0)

    non_periodic_mesh = lattice.generate_vtk(size=0.1, periodic=False)
    periodic_mesh = lattice.generate_vtk(size=0.1, periodic=True)

    assert non_periodic_mesh is not periodic_mesh
    assert not is_periodic(non_periodic_mesh.points)
    assert is_periodic(periodic_mesh.points)

