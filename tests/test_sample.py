"""Miscellaneous tests."""

import shutil
import warnings
from sys import platform

import microgen

USE_NEPER = shutil.which("neper") is not None
if not USE_NEPER:
    warnings.warn(
        "Neper is not installed. Please install Neper before running this command",
        stacklevel=2,
    )


def test_neper() -> None:
    """Test neper generation."""
    if USE_NEPER and platform != "win32":
        microgen.Neper.run(
            filename="tests/data/neper.tess",
            n_cells=2,
            cube_dim=(1, 1, 1),
        )
        microgen.Neper.voronoi_from_tess_file("tests/data/neper.tess")


def test_operations() -> None:
    """Test operations on shapes."""
    elem = microgen.Box(center=(0.5, 0.5, 0.5), dim=(1, 1, 1))
    shape1 = elem.generate_cad()
    phase1 = microgen.Phase.from_cad(shape1)

    elem = microgen.Box(center=(0, 0, 0), dim=(0.5, 0.5, 0.5))
    shape2 = elem.generate_cad()
    microgen.rescale(shape2, 2.0)

    microgen.cut_phase_by_shape_list(phase_to_cut=phase1, shapes=[shape2])

    microgen.cut_phases_by_shape(phases=[phase1], cut_obj=shape2)

    rve = microgen.Rve(dim=1, center=(0.5, 0.5, 0.5))
    microgen.repeat_shape(shape1, rve, grid=(2, 2, 2))
