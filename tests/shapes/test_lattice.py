"""Test the lattice shapes."""

from typing import Literal

import numpy as np
import pytest

pytest.importorskip("OCP")

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


# Each density fit runs ~10 root-finding steps, each a full CAD fusion of the
# lattice, so keep the grid small: the relative density does not depend on
# cell_size (checked on the cheap Cubic only), and OctetTruss (~1 min per fit)
# is covered by the fusion tests of test_cad_compat.py.
@pytest.mark.parametrize(
    ("shape", "expected_density", "cell_size"),
    [
        (Cubic, 0.25, 1.0),
        (Cubic, 0.8, 1.0),
        (Diamond, 0.25, 1.0),
        (Diamond, 0.8, 1.0),
        (Cubic, 0.25, 0.5),
        (Cubic, 0.25, 2.0),
    ],
)
def test_lattice_given_density_and_cell_size_must_match_computed_density(
    shape: Literal[Diamond, Cubic],
    expected_density: float,
    cell_size: float,
) -> None:
    """The radius fitted to ``density`` gives that density (CAD volume)."""
    # Arrange
    shape_fit_to_density = shape(density=expected_density, cell_size=cell_size)

    # Act
    shape_density = shape_fit_to_density.density

    # Assert
    assert np.isclose(expected_density, shape_density, rtol=0.01)


def test_lattice_mesh_for_fem_periodic_must_produce_periodic_mesh() -> None:
    """mesh_for_fem(periodic=True) must produce a mesh with periodic node positions."""
    lattice = OctetTruss(strut_radius=0.05, cell_size=1.0)
    mesh = lattice.mesh_for_fem(size=0.1, periodic=True)

    assert is_periodic(mesh.points), (
        "Mesh generated with periodic=True must be periodic"
    )


def test_lattice_mesh_for_fem_must_not_reuse_non_periodic_mesh_for_periodic_request() -> (
    None
):
    """mesh_for_fem must cache by parameters, not by instance only."""
    lattice = OctetTruss(strut_radius=0.05, cell_size=1.0)

    non_periodic_mesh = lattice.mesh_for_fem(size=0.1, periodic=False)
    periodic_mesh = lattice.mesh_for_fem(size=0.1, periodic=True)

    assert non_periodic_mesh is not periodic_mesh
    assert not is_periodic(non_periodic_mesh.points)
    assert is_periodic(periodic_mesh.points)
