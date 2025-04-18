"""Test the lattice shapes."""

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
