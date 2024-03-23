"""Test the periodic function."""

import numpy as np
import pytest

from microgen import Capsule, Phase, Rve, Sphere, periodic


def generate_sphere(x: float, y: float, z: float, rve: Rve):
    """Generate a sphere at the given coordinates."""
    elem = Sphere(center=(x, y, z), radius=0.1)
    phase = Phase(shape=elem.generate())
    periodic(phase=phase, rve=rve)


@pytest.mark.filterwarnings("ignore:Object intersecting")
def test_periodic():
    """Test the periodic function."""
    rve = Rve(dim=1, center=(0.5, 0.5, 0.5))

    # test x- and x+ faces intersected
    elem = Capsule(center=(0.5, 0, 0.5), height=1, radius=0.1)
    phase = Phase(shape=elem.generate())
    periodic(phase=phase, rve=rve)

    xyz = np.mgrid[3 * (slice(0.0, 1.0, 3),)]
    nodes = np.array([axis.flatten() for axis in xyz]).T
    for node in nodes:
        generate_sphere(*node, rve=rve)
