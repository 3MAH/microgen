"""Tests for Rve class."""

import numpy as np
import pytest

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


def test_rve_deprecated_dim() -> None:
    """Test Rve dim with deprecated inputs, dim_x, dim_y, dim_z are deprecated."""
    with pytest.warns(DeprecationWarning):
        rve = Rve(dim_x=1, dim_y=2, dim_z=1)
        assert np.all(rve.dim == [1, 2, 1])

    with pytest.warns(DeprecationWarning):
        rve = Rve(dim_x=2)
        assert np.all(rve.dim == [2, 1, 1])

    with pytest.warns(DeprecationWarning):
        rve = Rve(dim_x=1, dim_y=1, dim_z=1, dim=(2, 2, 2))
        assert np.all(rve.dim == [1, 1, 1])


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

    rve = Rve.from_min_max(x_min=-1, x_max=0, y_min=0, y_max=2, z_min=-3, z_max=3)
    assert np.all(rve.center == [-0.5, 1, 0])
    assert np.all(rve.dim == [1, 2, 6])
