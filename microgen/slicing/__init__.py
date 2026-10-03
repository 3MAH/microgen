"""Experimental meshless slicing, in millimeters and an explicit build frame.

Install ``microgen[slicing]`` for contour extraction and UVJ writing. PySLM is
optional and supplied by the caller. No printer-ready parameters are invented.
"""

from .engine import ImplicitContourSource, slice_contours, slice_rasters
from .field import BoundedSolid
from .pyslm import hatch_with_pyslm, to_pyslm_boundaries
from .types import (
    ContourLayer,
    ContourRegion,
    ContourSource,
    LayerPlan,
    LayerSpec,
    PixelGrid,
    RasterLayer,
)
from .uvj import read_uvj_profile, write_uvj

__all__ = [
    "BoundedSolid",
    "ContourLayer",
    "ContourRegion",
    "ContourSource",
    "ImplicitContourSource",
    "LayerPlan",
    "LayerSpec",
    "PixelGrid",
    "RasterLayer",
    "hatch_with_pyslm",
    "read_uvj_profile",
    "slice_contours",
    "slice_rasters",
    "to_pyslm_boundaries",
    "write_uvj",
]
