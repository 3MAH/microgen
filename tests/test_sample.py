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
            nbCell=2,
            dimCube=(1, 1, 1),
        )
        microgen.parseNeper("tests/data/neper")
        microgen.Neper.generateVoronoiFromTessFile("tests/data/neper.tess")


def test_operations() -> None:
    """Test operations on shapes."""
    elem = microgen.Box(center=(0.5, 0.5, 0.5), dim=(1, 1, 1))
    shape1 = elem.generate()
    phase1 = microgen.Phase(shape=shape1)

    elem = microgen.Box(center=(0, 0, 0), dim=(0.5, 0.5, 0.5))
    shape2 = elem.generate()
    microgen.rescale(shape2, 2.0)

    microgen.cutPhaseByShapeList(phaseToCut=phase1, cqShapeList=[shape2])

    microgen.cutPhasesByShape(phaseList=[phase1], cut_obj=shape2)

    rve = microgen.Rve(dim=1, center=(0.5, 0.5, 0.5))
    microgen.repeatShape(shape1, rve, grid=(2, 2, 2))
