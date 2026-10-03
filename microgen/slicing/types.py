"""Experimental manufacturing layer contracts, in millimeters and build coordinates."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]


def _finite(values: tuple[float, ...], name: str) -> None:
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"{name} must be finite")


@dataclass(frozen=True)
class LayerSpec:
    """One slab and its sampling plane; equal bounds describe a standalone plane."""

    index: int
    z_bottom: float
    z_top: float
    z_sample: float

    def __post_init__(self) -> None:
        _finite((self.z_bottom, self.z_top, self.z_sample), "Layer coordinates")
        if (
            not isinstance(self.index, int)
            or isinstance(self.index, bool)
            or self.index < 0
            or not self.z_bottom <= self.z_sample <= self.z_top
        ):
            raise ValueError("Invalid layer index or sampling plane")

    @property
    def thickness(self) -> float:
        """Slab thickness in millimeters."""
        return self.z_top - self.z_bottom

    @classmethod
    def plane(cls, z: float) -> LayerSpec:
        """Describe a single contour plane without assigning manufacturing thickness."""
        return cls(0, z, z, z)


@dataclass(frozen=True)
class LayerPlan:
    """Ordered, contiguous positive-thickness slabs; images are never stored here."""

    layers: tuple[LayerSpec, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "layers", tuple(self.layers))
        if not self.layers:
            raise ValueError("A layer plan must contain at least one layer")
        for index, layer in enumerate(self.layers):
            if layer.index != index or layer.thickness <= 0:
                raise ValueError(
                    "Layers must have consecutive indices and positive thickness"
                )
            if index and not math.isclose(
                self.layers[index - 1].z_top,
                layer.z_bottom,
                abs_tol=1e-10,
                rel_tol=0,
            ):
                raise ValueError("Layer slabs must be contiguous")

    @classmethod
    def uniform(cls, *, z_min: float, z_max: float, thickness: float) -> LayerPlan:
        """Use midpoint sampling, clipping the last slab to z_max."""
        _finite((z_min, z_max, thickness), "Layer plan parameters")
        if z_max <= z_min or thickness <= 0:
            raise ValueError("Require z_max > z_min and thickness > 0")
        count = math.ceil(math.nextafter((z_max - z_min) / thickness, -math.inf))
        layers = []
        for index in range(count):
            bottom = z_min + index * thickness
            top = min(z_max, z_min + (index + 1) * thickness)
            layers.append(LayerSpec(index, bottom, top, (bottom + top) / 2))
        return cls(tuple(layers))


@dataclass(frozen=True)
class PixelGrid:
    """Pixel extent; origin is lower left, image row zero is highest build Y.

    Pixel centers are half a pitch inside the extent. Columns increase along X,
    rows decrease along Y. Printer mirroring belongs to the printer adapter.
    """

    width_px: int
    height_px: int
    pitch_mm: tuple[float, float]
    origin_mm: tuple[float, float] = (0.0, 0.0)

    def __post_init__(self) -> None:
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value <= 0
            for value in (self.width_px, self.height_px)
        ):
            raise ValueError("Pixel dimensions must be positive integers")
        if len(self.pitch_mm) != 2 or len(self.origin_mm) != 2:
            raise ValueError("Pixel pitch and origin must each have two coordinates")
        _finite((*self.pitch_mm, *self.origin_mm), "Pixel grid parameters")
        if min(self.pitch_mm) <= 0:
            raise ValueError("Pixel pitches must be positive")
        object.__setattr__(self, "pitch_mm", tuple(self.pitch_mm))
        object.__setattr__(self, "origin_mm", tuple(self.origin_mm))

    @property
    def size_mm(self) -> tuple[float, float]:
        """Physical image width and height."""
        return self.width_px * self.pitch_mm[0], self.height_px * self.pitch_mm[1]


def _ring(points: FloatArray, *, clockwise: bool) -> FloatArray:
    ring = np.array(points, dtype=float, copy=True)
    if (
        ring.ndim != 2
        or ring.shape[1] != 2
        or len(ring) < 3
        or not np.isfinite(ring).all()
    ):
        raise ValueError("A ring requires at least three finite XY vertices")
    if not np.array_equal(ring[0], ring[-1]):
        ring = np.vstack((ring, ring[0]))
    # Translate before computing area to avoid cancellation at large origins.
    shifted = ring - ring[0]
    area = np.sum(shifted[:-1, 0] * shifted[1:, 1] - shifted[1:, 0] * shifted[:-1, 1])
    if area == 0:
        raise ValueError("A ring must enclose nonzero area")
    if (area < 0) != clockwise:
        ring = ring[::-1].copy()
    ring.setflags(write=False)
    return ring


@dataclass(frozen=True, eq=False)
class ContourRegion:
    """A solid polygon, with closed CCW exterior and closed CW holes.

    Construction validates polygon topology using the slicing extra. Arrays are
    copied so a caller cannot accidentally mutate a region through its inputs.
    """

    outer: FloatArray
    holes: tuple[FloatArray, ...] = ()
    part_id: int = 1

    def __post_init__(self) -> None:
        from shapely.geometry import Polygon

        object.__setattr__(self, "outer", _ring(self.outer, clockwise=False))
        object.__setattr__(
            self,
            "holes",
            tuple(_ring(hole, clockwise=True) for hole in self.holes),
        )
        polygon = Polygon(self.outer, self.holes)
        if (
            not polygon.is_valid
            or polygon.is_empty
            or not isinstance(self.part_id, int)
            or isinstance(self.part_id, bool)
            or self.part_id < 1
        ):
            raise ValueError("Invalid polygon topology or part ID")


@dataclass(frozen=True)
class ContourLayer:
    """Solid regions at a sampling plane, with explicit discretization diagnostics."""

    position: LayerSpec
    regions: tuple[ContourRegion, ...]
    sampling_pitch_mm: tuple[float, float] | None = None
    diagnostics: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "regions", tuple(self.regions))
        object.__setattr__(self, "diagnostics", tuple(self.diagnostics))
        if self.sampling_pitch_mm is not None:
            if len(self.sampling_pitch_mm) != 2:
                raise ValueError("Contour sampling pitch requires two coordinates")
            _finite(self.sampling_pitch_mm, "Contour sampling pitch")
            if min(self.sampling_pitch_mm) <= 0:
                raise ValueError("Contour sampling pitches must be positive")
            object.__setattr__(self, "sampling_pitch_mm", tuple(self.sampling_pitch_mm))


@dataclass(frozen=True, eq=False)
class RasterLayer:
    """One uint8 coverage image, with zero outside and 255 fully inside."""

    position: LayerSpec
    grid: PixelGrid
    image: npt.NDArray[np.uint8]
    supersampling: int = 1

    def __post_init__(self) -> None:
        image = np.asarray(self.image)
        if image.dtype != np.uint8 or image.shape != (
            self.grid.height_px,
            self.grid.width_px,
        ):
            raise ValueError(
                "Raster must be uint8 with the pixel grid's height and width"
            )
        if (
            not isinstance(self.supersampling, int)
            or isinstance(self.supersampling, bool)
            or self.supersampling < 1
        ):
            raise ValueError("Supersampling must be a positive integer")
        image = image.copy()
        image.setflags(write=False)
        object.__setattr__(self, "image", image)


class ContourSource(Protocol):
    """Produce geometry in millimeters in a common build frame.

    tolerance requests a spatial sampling scale, not certified topology or a
    guaranteed Hausdorff error. Implementations must report their limitations.
    """

    def contours_at(self, z: float, *, tolerance: float) -> ContourLayer:
        """Extract one plane, including any unresolved sampling diagnostics."""
        ...
