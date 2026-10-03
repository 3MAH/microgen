"""Optional real UVtools decoding/conversion test, enabled by UVTOOLS_CMD."""

import io
import json
import os
import subprocess
import zipfile

import numpy as np
import pytest

pytest.importorskip("PIL")
from PIL import Image

from microgen.slicing import (
    BoundedSolid,
    LayerPlan,
    PixelGrid,
    read_uvj_profile,
    slice_rasters,
    write_uvj,
)


def test_uvtools_native_goo_roundtrip(tmp_path):
    command = os.environ.get("UVTOOLS_CMD")
    if not command:
        pytest.skip("Set UVTOOLS_CMD to the UVtools command-line executable")
    grid = PixelGrid(32, 24, (0.1, 0.1))
    plan = LayerPlan.uniform(z_min=0, z_max=0.2, thickness=0.1)
    solid = BoundedSolid(
        lambda x, y, z: np.maximum(x - 1.5, 1 - y), (0, 3.2, 0, 2.4, 0, 0.2)
    )
    layers = list(slice_rasters(solid, plan, grid=grid, supersampling=2))
    # Test-only process values. UVtools normalizes this into a real profile fixture.
    profile = {
        "Properties": {
            "Size": {
                "X": 32,
                "Y": 24,
                "Millimeter": {"X": 3.2, "Y": 2.4},
                "Layers": 2,
                "LayerHeight": 0.1,
            },
            "Exposure": {
                "LightOnTime": 2,
                "LightPWM": 255,
                "LiftHeight": 5,
                "LiftSpeed": 60,
                "RetractSpeed": 120,
            },
            "Bottom": {
                "Count": 1,
                "LightOnTime": 20,
                "LightPWM": 255,
                "LiftHeight": 5,
                "LiftSpeed": 60,
                "RetractSpeed": 120,
            },
            "AntiAliasLevel": 1,
        }
    }
    seed = write_uvj(layers, tmp_path / "seed.uvj", profile=profile)

    def convert(source, target_type, destination):
        result = subprocess.run(
            [command, "convert", str(source), target_type, str(destination)],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=90,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert destination.is_file(), result.stdout + result.stderr

    normalized = tmp_path / "profile.uvj"
    convert(seed, "uvj", normalized)
    source = write_uvj(
        layers, tmp_path / "source.uvj", profile=read_uvj_profile(normalized)
    )
    native = tmp_path / "native.goo"
    convert(source, "goo", native)
    roundtrip = tmp_path / "roundtrip.uvj"
    convert(native, "uvj", roundtrip)
    with zipfile.ZipFile(roundtrip) as archive:
        metadata = json.loads(archive.read("config.json"))
        assert metadata["Properties"]["Size"]["X"] == grid.width_px
        assert metadata["Properties"]["Size"]["Y"] == grid.height_px
        assert metadata["Properties"]["Size"]["LayerHeight"] == pytest.approx(0.1)
        assert [layer["Z"] for layer in metadata["Layers"]] == pytest.approx([0.1, 0.2])
        assert [
            layer["Exposure"]["LightOnTime"] for layer in metadata["Layers"]
        ] == pytest.approx([20, 2])
        for index, layer in enumerate(layers):
            png = archive.read(f"slice/{index:08d}.png")
            np.testing.assert_array_equal(
                np.asarray(Image.open(io.BytesIO(png))), layer.image
            )
