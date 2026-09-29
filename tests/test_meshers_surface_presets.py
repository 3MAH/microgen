"""Direct surface presets use native meshers and preserve rigid placement."""

import meshers
import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from microgen import Tpms
from microgen.shape.surface_functions import gyroid

pytestmark = pytest.mark.skipif(
    not hasattr(meshers, "generate_surface")
    or not hasattr(meshers._meshers, "generate_surface"),
    reason="requires the experimental-surfaces meshers build",
)


@pytest.mark.parametrize("optimization", ["fast", "accurate", "quality"])
def test_surface_presets_match_native_generation(optimization):
    rotation = Rotation.from_euler("z", 31, degrees=True)
    shape = Tpms(
        gyroid, offset=0.6, resolution=12, center=(2, 3, 4), orientation=rotation
    )
    options = dict(refine_edges=optimization != "fast")
    if optimization != "quality":
        options.update(polish_passes=0, improvement_rounds=0, smoothing_iterations=0)
    expected = meshers.generate_surface(
        lambda x, y, z: shape.raw_field(x, y, z) / 0.3,
        bounds=shape._bounds,
        cells=11,
        band=(-1, 1),
        periodic=(True,) * 3,
        **options,
    )
    actual = shape.generate_meshers_surface(
        optimization=optimization, periodic=(True,) * 3
    )
    np.testing.assert_allclose(
        actual.points, rotation.apply(expected.points) + [2, 3, 4], atol=1e-12
    )
    np.testing.assert_array_equal(actual.triangles, expected.triangles)
    np.testing.assert_array_equal(actual.labels, expected.labels)
    assert actual.diagnostics["optimization"] == optimization


def test_surface_overrides_and_validation():
    shape = Tpms(gyroid, offset=0.6, resolution=12)
    a = shape.generate_meshers_surface(
        optimization="quality", polish_passes=0, improvement_rounds=0
    )
    b = shape.generate_meshers_surface(optimization="accurate")
    np.testing.assert_array_equal(a.points, b.points)
    np.testing.assert_array_equal(a.triangles, b.triangles)
    with pytest.raises(ValueError, match="optimization"):
        shape.generate_meshers_surface(optimization="unknown")
    with pytest.raises(NotImplementedError, match="equal grid"):
        Tpms(
            gyroid, offset=0.6, repeat_cell=(2, 1, 1), resolution=12
        ).generate_meshers_surface()
