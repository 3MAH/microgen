"""A polygon adapter for an existing PySLM hatcher, with no PySLM dependency."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Protocol, TypeVar

import numpy as np

from .types import ContourLayer, FloatArray, LayerSpec

HatchResult_co = TypeVar("HatchResult_co", covariant=True)


class Hatcher(Protocol[HatchResult_co]):
    """The PySLM operation consumed by this adapter."""

    def hatch(self, boundaryFeature: list[FloatArray]) -> HatchResult_co:
        """Use caller-configured scan settings to hatch closed polygon boundaries."""
        ...


def to_pyslm_boundaries(layer: ContourLayer) -> list[FloatArray]:
    """Return closed CCW exteriors and CW holes as independent writable arrays.

    PySLM consumes planar paths. Its native layer Z/index, model IDs, and
    build-style IDs remain the responsibility of the process adapter.
    """
    return [
        np.array(ring, copy=True)
        for region in layer.regions
        for ring in (region.outer, *region.holes)
    ]


def hatch_with_pyslm(
    layers: Iterable[ContourLayer],
    hatcher: Hatcher[HatchResult_co],
) -> Iterator[tuple[LayerSpec, HatchResult_co | None]]:
    """Stream positions and native hatcher results; preserve empty layers.

    Part IDs are not represented by PySLM's boundary-list argument. Reject
    mixed-part layers so their process identity cannot be silently flattened.
    """
    for layer in layers:
        if len({region.part_id for region in layer.regions}) > 1:
            raise ValueError("Hatch each part separately to retain process identity")
        boundaries = to_pyslm_boundaries(layer)
        yield layer.position, hatcher.hatch(boundaries) if boundaries else None
