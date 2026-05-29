"""Tests for Rve class."""

import dataclasses

import numpy as np
import pytest
import pyvista as pv

from microgen import Rve

# ruff: noqa: S101 assert https://docs.astral.sh/ruff/rules/assert/
# ruff: noqa: E501 line-too-long https://docs.astral.sh/ruff/rules/line-too-long/


def test_rve_center_must_have_expected_value() -> None:
    """Test Rve center with different input types."""
    rve = Rve()
    assert np.all(rve.center == [0, 0, 0])

    rve = Rve(center=(1, 2, 3))
    assert np.all(rve.center == [1, 2, 3])

    rve = Rve(center=[1, 2, 3])
    assert np.all(rve.center == [1, 2, 3])

    rve = Rve(center=np.array([1, 2, 3]))
    assert np.all(rve.center == [1, 2, 3])


def test_rve_invalid_center_value() -> None:
    """Test Rve center with invalid inputs, center must be an array or Sequence of length 3."""
    invalid_center_msg = "center must be an array or Sequence of length 3"
    with pytest.raises(ValueError, match=invalid_center_msg):
        Rve(center=(1, 2))

    with pytest.raises(ValueError, match=invalid_center_msg):
        Rve(center=[1, 2, 3, 4])

    with pytest.raises(ValueError, match=invalid_center_msg):
        Rve(center=np.array([1, 2]))  # type: ignore[arg-type]

    with pytest.raises(ValueError, match=invalid_center_msg):
        Rve(center=1)  # type: ignore[arg-type]


def test_rve_dim_must_have_expected_value() -> None:
    """Test Rve dim with different input types."""
    rve = Rve()
    assert np.all(rve.dim == [1, 1, 1])

    rve = Rve(dim=2)
    assert np.all(rve.dim == [2, 2, 2])

    rve = Rve(dim=(2, 3, 4))
    assert np.all(rve.dim == [2, 3, 4])

    rve = Rve(dim=[2, 3, 4])
    assert np.all(rve.dim == [2, 3, 4])

    rve = Rve(dim=np.array([2, 3, 4]))
    assert np.all(rve.dim == [2, 3, 4])


def test_rve_given_negative_or_zero_dim_must_raise_error() -> None:
    """Test Rve dim with negative values must raise error."""
    invalid_dim_msg = "dimensions of the RVE must be greater than 0"
    with pytest.raises(ValueError, match=invalid_dim_msg):
        Rve(dim=-1)

    with pytest.raises(ValueError, match=invalid_dim_msg):
        Rve(dim=0)

    with pytest.raises(ValueError, match=invalid_dim_msg):
        Rve(dim=(1, -1, 0))


def test_rve_from_min_max_must_return_expected_rve() -> None:
    """Test Rve from_min_max class method must return expected Rve."""
    rve = Rve.from_min_max()
    assert np.all(rve.center == [0, 0, 0])
    assert np.all(rve.dim == [1, 1, 1])

    rve = Rve.from_min_max(min_point=(-1, 0, -3), max_point=(0, 2, 3))
    assert np.all(rve.center == [-0.5, 1, 0])
    assert np.all(rve.dim == [1, 2, 6])


def test_rve_pbc_default_is_fully_periodic() -> None:
    """``pbc`` defaults to ``(True, True, True)``."""
    assert Rve().pbc == (True, True, True)


def test_rve_pbc_scalar_broadcasts_to_all_axes() -> None:
    """A single bool is broadcast to all three axes."""
    assert Rve(pbc=False).pbc == (False, False, False)
    assert Rve(pbc=True).pbc == (True, True, True)


def test_rve_pbc_accepts_per_axis_tuple() -> None:
    """A length-3 sequence sets each axis independently."""
    assert Rve(pbc=(True, False, True)).pbc == (True, False, True)
    assert Rve(pbc=[False, False, True]).pbc == (False, False, True)


def test_rve_invalid_pbc_raises() -> None:
    """``pbc`` of wrong length or type raises ``ValueError``."""
    with pytest.raises(ValueError, match="pbc must be a bool"):
        Rve(pbc=(True, False))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="pbc must be a bool"):
        Rve(pbc=1.5)  # type: ignore[arg-type]


def test_rve_is_frozen() -> None:
    """Reassigning ``center``/``dim``/``pbc`` after construction raises FrozenInstanceError."""
    rve = Rve(dim=2)
    with pytest.raises(dataclasses.FrozenInstanceError):
        rve.center = np.array([1.0, 2.0, 3.0])  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        rve.dim = np.array([3.0, 3.0, 3.0])  # type: ignore[misc]
    with pytest.raises(dataclasses.FrozenInstanceError):
        rve.pbc = (False, False, False)  # type: ignore[misc]


def test_rve_min_max_points_are_consistent() -> None:
    """``min_point`` / ``max_point`` are derived from center ± 0.5·dim."""
    rve = Rve(center=(1.0, 2.0, 3.0), dim=(2.0, 4.0, 6.0))
    assert np.allclose(rve.min_point, [0.0, 0.0, 0.0])
    assert np.allclose(rve.max_point, [2.0, 4.0, 6.0])


def test_rve_grid_with_scalar_resolution() -> None:
    """``Rve.grid(n)`` returns a StructuredGrid spanning the cell with n^3 points."""
    rve = Rve(center=(0.0, 0.0, 0.0), dim=(2.0, 4.0, 6.0))
    grid = rve.grid(5)
    assert isinstance(grid, pv.StructuredGrid)
    assert grid.dimensions == (5, 5, 5)
    pts = np.asarray(grid.points)
    assert np.isclose(pts[:, 0].min(), -1.0)
    assert np.isclose(pts[:, 0].max(), 1.0)
    assert np.isclose(pts[:, 1].min(), -2.0)
    assert np.isclose(pts[:, 1].max(), 2.0)
    assert np.isclose(pts[:, 2].min(), -3.0)
    assert np.isclose(pts[:, 2].max(), 3.0)


def test_rve_grid_with_per_axis_resolution() -> None:
    """``Rve.grid((nx, ny, nz))`` honors per-axis resolution."""
    rve = Rve(center=(0.0, 0.0, 0.0), dim=1.0)
    grid = rve.grid((3, 4, 5))
    assert grid.dimensions == (3, 4, 5)


def test_rve_grid_invalid_resolution_raises() -> None:
    """Wrong-length resolution raises ``ValueError``."""
    rve = Rve()
    with pytest.raises(ValueError, match="resolution must be an int"):
        rve.grid((3, 4))  # type: ignore[arg-type]


def test_rve_repr_round_trips_inputs() -> None:
    """``repr(rve)`` contains the canonical ``center``/``dim``/``pbc`` triple."""
    rve = Rve(center=(1.0, 2.0, 3.0), dim=(2.0, 2.0, 2.0), pbc=(True, False, True))
    text = repr(rve)
    assert "center=(1.0, 2.0, 3.0)" in text
    assert "dim=(2.0, 2.0, 2.0)" in text
    assert "pbc=(True, False, True)" in text
