"""Analytic geometry, streaming, and receiver contracts for experimental slicing."""

import io
import json
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("skimage")
pytest.importorskip("shapely")
pytest.importorskip("PIL")

from PIL import Image
from shapely.geometry import Point, Polygon

from microgen.shape.implicit_ops import from_field
from microgen.slicing import (
    BoundedSolid,
    ContourLayer,
    ContourRegion,
    ImplicitContourSource,
    LayerPlan,
    LayerSpec,
    PixelGrid,
    RasterLayer,
    hatch_with_pyslm,
    read_uvj_profile,
    slice_contours,
    slice_rasters,
    to_pyslm_boundaries,
    write_uvj,
)


def sphere(radius=1):
    return BoundedSolid(
        lambda x, y, z: x * x + y * y + z * z - radius**2, (-2, 2, -2, 2, -2, 2)
    )


def area(layer):
    return sum(Polygon(region.outer, region.holes).area for region in layer.regions)


@pytest.fixture
def profile():
    # Explicit test process settings, not recommended printer parameters.
    return {
        "Properties": {
            "Size": {
                "X": 4,
                "Y": 4,
                "Millimeter": {"X": 4.0, "Y": 4.0},
                "Layers": 100,
                "LayerHeight": 0.1,
            },
            "Exposure": {"LightOnTime": 2.0, "LiftHeight": 5.0, "LightPWM": 255},
            "Bottom": {
                "Count": 1,
                "LightOnTime": 20.0,
                "LiftHeight": 6.0,
                "LightPWM": 255,
            },
            "AntiAliasLevel": 1,
            "Vendor": {"test": {"retained": True}},
        },
        "Layers": [{"Z": 99, "Exposure": {"LightOnTime": 99}}],
    }


def rasters():
    plan = LayerPlan.uniform(z_min=0, z_max=0.2, thickness=0.1)
    grid = PixelGrid(4, 4, (1.0, 1.0))
    solid = BoundedSolid(lambda x, y, z: x + y - 3, (0, 4, 0, 4, 0, 0.2))
    return slice_rasters(solid, plan, grid=grid)


def test_layer_plan_midpoints_and_last_partial_slab():
    plan = LayerPlan.uniform(z_min=2, z_max=2.25, thickness=0.1)
    assert [p.index for p in plan.layers] == [0, 1, 2]
    assert [p.z_sample for p in plan.layers] == pytest.approx([2.05, 2.15, 2.225])
    assert plan.layers[-1].thickness == pytest.approx(0.05)
    assert len(LayerPlan.uniform(z_min=0, z_max=0.3, thickness=0.1).layers) == 3


@pytest.mark.parametrize(
    "arguments",
    [
        {"z_min": 1, "z_max": 0, "thickness": 0.1},
        {"z_min": 0, "z_max": 1, "thickness": 0},
        {"z_min": 0, "z_max": float("nan"), "thickness": 0.1},
    ],
)
def test_invalid_layer_plans(arguments):
    with pytest.raises(ValueError):
        LayerPlan.uniform(**arguments)


def test_plan_rejects_gaps_and_reordered_layers():
    with pytest.raises(ValueError):
        LayerPlan((LayerSpec(0, 0, 1, 0.5), LayerSpec(1, 2, 3, 2.5)))
    with pytest.raises(ValueError):
        LayerPlan((LayerSpec(1, 0, 1, 0.5),))


def test_field_only_shape_adapter_never_meshes(monkeypatch):
    shape = from_field(
        lambda x, y, z: x * x + y * y + z * z - 1, bounds=(-2, 2, -2, 2, -2, 2)
    )

    def forbid(*args, **kwargs):
        raise AssertionError("A mesh or 3D grid was requested")

    monkeypatch.setattr(shape, "generate_surface_mesh", forbid)
    monkeypatch.setattr(shape, "_sample_implicit_grid", forbid)
    solid = BoundedSolid.from_shape(shape)
    layer = ImplicitContourSource(solid).contours_at(0, tolerance=0.05)
    assert area(layer) == pytest.approx(np.pi, abs=0.02)


def test_transform_and_bounded_field():
    transform = np.eye(4)
    transform[:3, 3] = (10, 20, 30)
    solid = BoundedSolid(
        lambda x, y, z: -np.ones_like(x), (-1, 1, -2, 2, -3, 3), transform
    )
    assert solid.build_bounds == (9, 11, 18, 22, 27, 33)
    assert solid.evaluate(
        np.array([10.0, 12.0]), np.array([20.0]), np.array([30.0])
    ).tolist() == [-1, 1]
    transform[:3, 3] = 0
    assert solid.build_bounds == (9, 11, 18, 22, 27, 33)


def test_rotated_solid_contours_use_build_coordinates():
    transform = np.array([[0.0, -1, 0, 5], [1, 0, 0, 6], [0, 0, 1, 0], [0, 0, 0, 1]])
    solid = BoundedSolid(
        lambda x, y, z: -np.ones_like(x), (-1, 1, -2, 2, -1, 1), transform
    )
    layer = ImplicitContourSource(solid).contours_at(0, tolerance=0.1)
    polygon = Polygon(layer.regions[0].outer)
    assert polygon.area == pytest.approx(8)
    assert polygon.bounds == pytest.approx((3, 5, 7, 7))


def test_nonfinite_fields_and_singular_transform_rejected():
    with pytest.raises(ValueError, match="invertible"):
        BoundedSolid(
            lambda x, y, z: -1,
            (-1, 1, -1, 1, -1, 1),
            np.zeros((4, 4)) + np.diag([0, 0, 0, 1]),
        )
    solid = BoundedSolid(lambda x, y, z: np.nan, (-1, 1, -1, 1, -1, 1))
    with pytest.raises(ValueError, match="nonfinite"):
        solid.evaluate(np.array([0]), np.array([0]), np.array([0]))


def test_sphere_section_matches_analytic_radius():
    z = 0.6
    layer = ImplicitContourSource(sphere()).contours_at(z, tolerance=0.025)
    assert len(layer.regions) == 1
    radii = np.linalg.norm(layer.regions[0].outer, axis=1)
    assert radii == pytest.approx(np.full_like(radii, 0.8), abs=0.0002)
    assert area(layer) == pytest.approx(np.pi * 0.8**2, abs=0.003)
    assert layer.diagnostics


def test_annulus_holes_and_nested_island():
    # Solid annulus plus a separate central disk.
    field = lambda x, y, z: np.minimum(
        np.maximum(np.hypot(x, y) - 1, 0.6 - np.hypot(x, y)), np.hypot(x, y) - 0.2
    )
    solid = BoundedSolid(field, (-1.2, 1.2, -1.2, 1.2, -1, 1))
    layer = ImplicitContourSource(solid).contours_at(0, tolerance=0.025)
    assert len(layer.regions) == 2
    assert sum(len(region.holes) for region in layer.regions) == 1
    assert area(layer) == pytest.approx(np.pi * (1 - 0.6**2 + 0.2**2), abs=0.005)


def test_full_box_and_empty_plane():
    solid = BoundedSolid(lambda x, y, z: -1, (0, 2, 0, 3, -1, 1))
    source = ImplicitContourSource(solid)
    assert area(source.contours_at(0, tolerance=0.2)) == pytest.approx(6)
    assert source.contours_at(2, tolerance=0.2).regions == ()
    assert source.contours_at(1, tolerance=0.2).regions == ()


def test_contour_allocation_guard():
    with pytest.raises(ValueError, match="max_grid_points"):
        ImplicitContourSource(sphere(), max_grid_points=100).contours_at(
            0, tolerance=0.001
        )


def test_raster_pixel_centers_and_image_orientation():
    grid = PixelGrid(3, 2, (1.0, 2.0), (10.0, 20.0))
    solid = BoundedSolid(
        lambda x, y, z: np.maximum(x - 11, 22 - y), (10, 13, 20, 24, -1, 1)
    )
    layer = next(
        slice_rasters(
            solid, LayerPlan.uniform(z_min=0, z_max=0.1, thickness=0.1), grid=grid
        )
    )
    np.testing.assert_array_equal(layer.image, [[255, 0, 0], [0, 0, 0]])


def test_supersampling_averages_occupancy():
    solid = BoundedSolid(lambda x, y, z: x - 0.5, (0, 1, 0, 1, -1, 1))
    plan = LayerPlan.uniform(z_min=0, z_max=0.1, thickness=0.1)
    layer = next(
        slice_rasters(solid, plan, grid=PixelGrid(1, 1, (1.0, 1.0)), supersampling=2)
    )
    assert layer.image[0, 0] == 128


def test_raster_evaluation_streams_tiles_and_layers():
    calls = []

    def field(x, y, z):
        calls.append(len(x))
        return -np.ones_like(x)

    solid = BoundedSolid(field, (0, 5, 0, 4, 0, 1))
    stack = slice_rasters(
        solid,
        LayerPlan.uniform(z_min=0, z_max=1, thickness=0.5),
        grid=PixelGrid(5, 4, (1.0, 1.0)),
        tile_rows=1,
    )
    assert calls == []
    first = next(stack)
    assert len(calls) == 4 and max(calls) == 5
    second = next(stack)
    assert len(calls) == 8
    assert not np.shares_memory(first.image, second.image)
    assert not first.image.flags.writeable


def test_contour_source_interface_streams_and_preserves_positions():
    calls = []

    class Source:
        def contours_at(self, z, *, tolerance):
            calls.append((z, tolerance))
            return ContourLayer(LayerSpec.plane(z), ())

    plan = LayerPlan.uniform(z_min=0, z_max=0.2, thickness=0.1)
    layers = slice_contours(Source(), plan, tolerance=0.01)
    assert calls == []
    assert next(layers).position == plan.layers[0]
    assert calls == [(0.05, 0.01)]


def test_region_owns_arrays_and_validates_topology():
    points = np.array([[0.0, 0], [2, 0], [2, 2], [0, 2]])
    region = ContourRegion(points)
    points[:] = 9
    assert region.outer[0].tolist() == [0, 0]
    assert np.array_equal(region.outer[0], region.outer[-1])
    assert not region.outer.flags.writeable
    with pytest.raises(ValueError):
        ContourRegion(
            np.array([[0.0, 0], [2, 0], [2, 2], [0, 2]]),
            (np.array([[3.0, 3], [4, 3], [4, 4]]),),
        )


def test_uvj_roundtrip_pixels_positions_profile_and_bottom_settings(tmp_path, profile):
    expected = list(rasters())
    destination = write_uvj(iter(expected), tmp_path / "test.uvj", profile=profile)
    with zipfile.ZipFile(destination) as archive:
        config = json.loads(archive.read("config.json"))
        assert config["Properties"]["Size"]["Layers"] == 2
        assert config["Properties"]["Vendor"]["test"]["retained"]
        assert [item["Z"] for item in config["Layers"]] == pytest.approx([0.1, 0.2])
        assert [item["Exposure"]["LightOnTime"] for item in config["Layers"]] == [20, 2]
        for index, layer in enumerate(expected):
            decoded = np.asarray(
                Image.open(io.BytesIO(archive.read(f"slice/{index:08d}.png")))
            )
            np.testing.assert_array_equal(decoded, layer.image)
    assert profile["Properties"]["Size"]["Layers"] == 100
    assert "Layers" not in read_uvj_profile(destination)


def test_uvj_rejects_wrong_printer_and_preserves_existing_file(tmp_path, profile):
    destination = tmp_path / "test.uvj"
    destination.write_bytes(b"previous")
    profile["Properties"]["Size"]["X"] = 9
    with pytest.raises(ValueError, match="resolution"):
        write_uvj(rasters(), destination, profile=profile)
    assert destination.read_bytes() == b"previous"
    assert list(tmp_path.iterdir()) == [destination]


def test_uvj_rejects_partial_last_layer_and_empty_stack(tmp_path, profile):
    plan = LayerPlan.uniform(z_min=0, z_max=0.15, thickness=0.1)
    solid = BoundedSolid(lambda x, y, z: -1, (0, 4, 0, 4, 0, 1))
    layers = slice_rasters(solid, plan, grid=PixelGrid(4, 4, (1.0, 1.0)))
    with pytest.raises(ValueError, match="constant layer height"):
        write_uvj(layers, tmp_path / "partial.uvj", profile=profile)
    with pytest.raises(ValueError, match="empty"):
        write_uvj([], tmp_path / "empty.uvj", profile=profile)
    assert not list(tmp_path.iterdir())


def test_real_pyslm_hatcher_respects_annulus_hole():
    pyslm = pytest.importorskip("pyslm.hatching")
    solid = BoundedSolid(
        lambda x, y, z: np.maximum(np.hypot(x, y) - 1, 0.4 - np.hypot(x, y)),
        (-1.2, 1.2, -1.2, 1.2, -1, 1),
    )
    layer = ImplicitContourSource(solid).contours_at(0, tolerance=0.025)
    hatcher = pyslm.Hatcher()
    hatcher.hatchDistance = 0.1
    hatcher.volumeOffsetHatch = 0.02
    hatcher.spotCompensation = 0.01
    position, result = next(hatch_with_pyslm([layer], hatcher))
    assert position == layer.position
    assert result is not None
    hatches = result.getHatchGeometry()
    assert hatches
    domain = Polygon(layer.regions[0].outer, layer.regions[0].holes).buffer(1e-5)
    from shapely.geometry import LineString

    for geometry in hatches:
        for segment in geometry.coords.reshape(-1, 2, 2):
            assert domain.covers(LineString(segment))
            assert not Point(0, 0).buffer(0.35).intersects(LineString(segment))
    assert all(ring.flags.writeable for ring in to_pyslm_boundaries(layer))


def test_hatching_retains_empty_layers_without_calling_backend():
    class NoHatching:
        def hatch(self, boundaryFeature):
            raise AssertionError("Empty layers must not be hatched")

    layer = ContourLayer(LayerSpec(0, 0, 0.1, 0.05), ())
    assert list(hatch_with_pyslm([layer], NoHatching())) == [(layer.position, None)]


def test_gyroid_example_uses_field_only_path(tmp_path, profile):
    seed = write_uvj(rasters(), tmp_path / "profile.uvj", profile=profile)
    output = tmp_path / "gyroid.uvj"
    example = (
        Path(__file__).resolve().parents[1]
        / "examples/ImplicitSlicing/gyroid_layers.py"
    )
    result = subprocess.run(
        [
            sys.executable,
            str(example),
            "--extent",
            "2",
            "--cell-size",
            "1",
            "--profile",
            str(seed),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "no mesh or 3D grid constructed" in result.stdout
    with zipfile.ZipFile(output) as archive:
        assert (
            json.loads(archive.read("config.json"))["Properties"]["Size"]["Layers"]
            == 20
        )
