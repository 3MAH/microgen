"""Tests for the grading of TPMS shapes."""

from __future__ import annotations

import numpy as np
import pytest
import pyvista as pv

from microgen import NormedDistance


@pytest.mark.parametrize(
    ("boundary_offset", "furthest_offset", "boundary_weight"),
    [
        (0.0, 3.0, 1.0),
        (3.0, 0.0, 0.5),
        (0.0, 3.0, 2.0),
    ],
)
def test_normed_distance_offsets_must_correspond(
    boundary_offset: float,
    furthest_offset: float,
    boundary_weight: float,
) -> None:
    """Test the NormedDistance grading function."""
    sphere = pv.Sphere()

    xrng = np.linspace(sphere.bounds[0], sphere.bounds[1], 5, endpoint=True)
    yrng = np.linspace(sphere.bounds[2], sphere.bounds[3], 5, endpoint=True)
    zrng = np.linspace(sphere.bounds[4], sphere.bounds[5], 5, endpoint=True)
    x, y, z = np.meshgrid(xrng, yrng, zrng)
    grid = pv.StructuredGrid(x, y, z)

    grading = NormedDistance(
        obj=sphere,
        boundary_offset=boundary_offset,
        furthest_offset=furthest_offset,
        boundary_weight=boundary_weight,
    )

    offset = grading.compute_offset(grid)

    middle = len(offset) // 2
    assert np.isclose(offset[middle], furthest_offset)
    assert np.isclose(offset[0], boundary_offset)
    assert np.isclose(offset[-1], boundary_offset)
