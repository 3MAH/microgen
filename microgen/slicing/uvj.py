"""Stream raster layers into UVtools' existing UVJ format using an explicit profile."""

from __future__ import annotations

import io
import json
import math
import os
import tempfile
import zipfile
from collections.abc import Iterable, Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

from .types import PixelGrid, RasterLayer


def read_uvj_profile(path: str | Path) -> dict[str, Any]:
    """Read printer/process metadata from a UVJ exported by UVtools.

    Global normal/bottom settings are retained. Source images and source
    per-layer overrides are not transferred to the new geometry.
    """
    with zipfile.ZipFile(path) as archive:
        config = json.loads(archive.read("config.json"))
    if not isinstance(config, dict) or not isinstance(config.get("Properties"), dict):
        raise TypeError("UVJ profile requires a Properties object")
    config.pop("Layers", None)
    return config


def _check_profile(config: dict[str, Any], layer: RasterLayer) -> None:
    try:
        properties = config["Properties"]
        size = properties["Size"]
        resolution = (size["X"], size["Y"])
        physical = (size["Millimeter"]["X"], size["Millimeter"]["Y"])
        height = float(size["LayerHeight"])
        exposure = properties["Exposure"]
        bottom = properties["Bottom"]
        normal_time = float(exposure["LightOnTime"])
        bottom_time = float(bottom["LightOnTime"])
        bottom_count = bottom["Count"]
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            "Require a UVtools UVJ printer profile with size and exposure settings"
        ) from error
    grid = layer.grid
    if resolution != (grid.width_px, grid.height_px) or max(resolution) > 65535:
        raise ValueError("Raster resolution must match the UVJ printer profile")
    if not all(
        math.isclose(float(a), b, abs_tol=1e-6, rel_tol=0)
        for a, b in zip(physical, grid.size_mm)
    ):
        raise ValueError("Physical raster size must match the UVJ printer profile")
    if grid.origin_mm != (0.0, 0.0):
        raise ValueError("UVJ requires a full-display raster with origin (0, 0)")
    if not math.isfinite(height) or height <= 0:
        raise ValueError("Profile layer height must be finite and positive")
    if (
        not isinstance(bottom_count, int)
        or isinstance(bottom_count, bool)
        or bottom_count < 0
    ):
        raise ValueError("Profile bottom-layer count must be a nonnegative integer")
    if not math.isfinite(normal_time) or normal_time <= 0:
        raise ValueError(
            "Require an explicit positive normal exposure from the printer profile"
        )
    if bottom_count and (not math.isfinite(bottom_time) or bottom_time <= 0):
        raise ValueError("Require positive bottom exposure for bottom layers")


def write_uvj(
    layers: Iterable[RasterLayer],
    path: str | Path,
    *,
    profile: Mapping[str, Any],
) -> Path:
    """Write PNG layers and UVJ metadata without retaining the image stack.

    The profile must match the complete display resolution, physical dimensions,
    and constant layer height. Slabs start at Z=0. Global exposure settings are
    applied by normal/bottom layer count, never inferred from geometry. Printer
    acceptance and exposure calibration remain external validation steps.

    Output replaces the destination only after the whole stack validates.
    """
    from PIL import Image

    destination = Path(path)
    config = deepcopy(dict(profile))
    positions: list[dict[str, Any]] = []
    first_grid: PixelGrid | None = None
    previous_top = 0.0
    with tempfile.NamedTemporaryFile(
        dir=destination.parent, suffix=".uvj", delete=False
    ) as handle:
        temporary = Path(handle.name)
    try:
        with zipfile.ZipFile(
            temporary, "w", compression=zipfile.ZIP_DEFLATED
        ) as archive:
            for index, layer in enumerate(layers):
                if first_grid is None:
                    _check_profile(config, layer)
                    first_grid = layer.grid
                properties = config["Properties"]
                height = float(properties["Size"]["LayerHeight"])
                if layer.grid != first_grid:
                    raise ValueError("All raster layers must use the same pixel grid")
                position = layer.position
                if position.index != index or not math.isclose(
                    position.z_bottom,
                    previous_top,
                    abs_tol=1e-7,
                    rel_tol=0,
                ):
                    raise ValueError(
                        "UVJ slabs must start at zero and be contiguous and ordered"
                    )
                if not math.isclose(
                    position.thickness, height, abs_tol=1e-7, rel_tol=0
                ):
                    raise ValueError(
                        "UVJ layers must match the profile's constant layer height"
                    )
                exposure = (
                    properties["Bottom"]
                    if index < properties["Bottom"]["Count"]
                    else properties["Exposure"]
                )
                exposure = {
                    key: value for key, value in exposure.items() if key != "Count"
                }
                positions.append({"Z": position.z_top, "Exposure": exposure})
                previous_top = position.z_top
                buffer = io.BytesIO()
                Image.fromarray(layer.image).save(buffer, format="PNG")
                archive.writestr(f"slice/{index:08d}.png", buffer.getvalue())
            if first_grid is None:
                raise ValueError("Cannot write an empty raster stack")
            config["Properties"]["Size"]["Layers"] = len(positions)
            config["Layers"] = positions
            archive.writestr(
                "config.json", json.dumps(config, indent=2, allow_nan=False)
            )
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination
