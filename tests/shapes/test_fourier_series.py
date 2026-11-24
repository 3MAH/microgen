"""Tests for the Fourier-series surface function."""

from __future__ import annotations

import numpy as np
import pytest

from microgen.shape import Tpms, surface_functions
from microgen.shape.surface_functions import fourier_series

GYROID_MODES = [(1, 1, 0), (1, -1, 0), (0, 1, 1), (0, 1, -1), (1, 0, 1), (-1, 0, 1)]


def _gyroid_series():
    return fourier_series(
        modes=GYROID_MODES,
        amplitudes=[0.5] * 6,
        phases=[-0.5 * np.pi] * 6,
    )


def _random_points(n: int = 500) -> np.ndarray:
    return np.random.default_rng(0).uniform(-5.0, 5.0, (3, n))


def test_fourier_series_matches_closed_form_gyroid():
    """The six Fourier modes of the gyroid reproduce the closed form."""
    x, y, z = _random_points()
    np.testing.assert_allclose(
        _gyroid_series()(x, y, z),
        surface_functions.gyroid(x, y, z),
        atol=1e-12,
    )


def test_fourier_series_default_phases_and_constant():
    """Phases default to zero and the constant term is added to the series."""
    x, y, z = _random_points()
    field = fourier_series(
        modes=[(1, 0, 0), (0, 2, 0), (1, 1, -1)],
        amplitudes=[1.0, -0.5, 0.25],
        constant=0.3,
    )
    expected = np.cos(x) - 0.5 * np.cos(2 * y) + 0.25 * np.cos(x + y - z) + 0.3
    np.testing.assert_allclose(field(x, y, z), expected, atol=1e-12)


def test_fourier_series_is_periodic():
    """Integer modes make the field 2*pi-periodic along each axis."""
    rng = np.random.default_rng(1)
    modes = rng.integers(-3, 4, size=(12, 3))
    field = fourier_series(
        modes=modes,
        amplitudes=rng.normal(size=12),
        phases=rng.uniform(0.0, 2.0 * np.pi, size=12),
    )
    x, y, z = _random_points()
    reference = field(x, y, z)
    for shift in np.eye(3) * 2.0 * np.pi:
        np.testing.assert_allclose(
            field(x + shift[0], y + shift[1], z + shift[2]),
            reference,
            atol=1e-10,
        )


def test_tpms_from_fourier_series_matches_gyroid():
    """A Tpms built on the series gives the same sheet and F-rep field."""
    series = Tpms(surface_function=_gyroid_series(), offset=0.3, resolution=20)
    closed = Tpms(surface_function=surface_functions.gyroid, offset=0.3, resolution=20)

    assert series.sheet.volume == pytest.approx(closed.sheet.volume, rel=1e-10)

    # The series is not autograd-traceable: the SDF normalization falls back
    # to finite differences, hence the looser tolerance.
    x, y, z = 0.05 * _random_points()
    np.testing.assert_allclose(
        series.as_sheet().field(x, y, z),
        closed.as_sheet().field(x, y, z),
        atol=1e-6,
    )


def test_tpms_from_fourier_series_fits_density():
    """The density solver works on a generic, non-minimal periodic surface."""
    field = fourier_series(
        modes=[(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (2, 1, 0)],
        amplitudes=[1.0, 1.0, 1.0, 0.8, 0.5],
        phases=[0.0, 0.0, 0.0, 1.0, -1.5],
    )
    geometry = Tpms(surface_function=field, density=0.3, resolution=30)
    assert geometry.sheet.volume == pytest.approx(0.3, rel=1e-2)


@pytest.mark.parametrize(
    ("modes", "amplitudes", "phases", "match"),
    [
        ([(1, 0)], [1.0], None, "shape"),
        ([(1, 0, 0.5)], [1.0], None, "integers"),
        ([(1, 0, 0), (0, 1, 0)], [1.0], None, "amplitudes"),
        ([(1, 0, 0)], [1.0], [0.0, 0.0], "phases"),
    ],
)
def test_fourier_series_rejects_invalid_input(modes, amplitudes, phases, match):
    """Malformed modes, amplitudes or phases raise ValueError."""
    with pytest.raises(ValueError, match=match):
        fourier_series(modes=modes, amplitudes=amplitudes, phases=phases)
