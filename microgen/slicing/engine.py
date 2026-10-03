"""Streaming raster sampling and a bounded uniform contouring baseline."""

from __future__ import annotations

import math
from collections.abc import Iterator
from dataclasses import dataclass, replace

import numpy as np

from .field import BoundedSolid
from .types import (
    ContourLayer,
    ContourRegion,
    ContourSource,
    FloatArray,
    LayerPlan,
    LayerSpec,
    PixelGrid,
    RasterLayer,
)


def slice_contours(
    source: ContourSource,
    plan: LayerPlan,
    *,
    tolerance: float,
) -> Iterator[ContourLayer]:
    """Attach slab positions to each source result, evaluating one layer at a time."""
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Contour tolerance must be finite and positive")
    for position in plan.layers:
        layer = source.contours_at(position.z_sample, tolerance=tolerance)
        if not math.isclose(
            layer.position.z_sample, position.z_sample, abs_tol=1e-10, rel_tol=0
        ):
            raise ValueError("Contour source returned the wrong sampling plane")
        yield replace(layer, position=position)


def slice_rasters(
    solid: BoundedSolid,
    plan: LayerPlan,
    *,
    grid: PixelGrid,
    supersampling: int = 1,
    tile_rows: int = 128,
) -> Iterator[RasterLayer]:
    """Sample subpixel centers and average occupancy, never scalar field values.

    This estimates planar area coverage, not exposure response or slab occupancy.
    Working field memory is limited to one tile; output memory is one image.
    """
    for name, value in (("Supersampling", supersampling), ("Tile rows", tile_rows)):
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    px, py = grid.pitch_mm
    ox, oy = grid.origin_mm
    cols = np.arange(grid.width_px, dtype=np.float64)[None, :]
    for position in plan.layers:
        image = np.empty((grid.height_px, grid.width_px), dtype=np.uint8)
        for start in range(0, grid.height_px, tile_rows):
            stop = min(start + tile_rows, grid.height_px)
            rows = np.arange(start, stop, dtype=np.float64)[:, None]
            coverage = np.zeros((stop - start, grid.width_px), dtype=float)
            for sy in range(supersampling):
                y = oy + (grid.height_px - rows - (sy + 0.5) / supersampling) * py
                for sx in range(supersampling):
                    x = ox + (cols + (sx + 0.5) / supersampling) * px
                    coverage += solid.evaluate(x, y, np.asarray(position.z_sample)) < 0
            image[start:stop] = np.rint(255 * coverage / supersampling**2).astype(
                np.uint8
            )
        yield RasterLayer(position, grid, image, supersampling)


def _refine_roots(
    solid: BoundedSolid,
    xy: FloatArray,
    *,
    pitch: tuple[float, float],
    z: float,
    origin: tuple[float, float],
) -> FloatArray:
    """Refine marching-square edge crossings with vectorized sign brackets."""
    cell = (xy - origin) / pitch
    vertical = np.isclose(cell[:, 0], np.rint(cell[:, 0]), atol=1e-8, rtol=0)
    low = cell.copy()
    high = cell.copy()
    low[:, 0] = np.where(vertical, np.rint(cell[:, 0]), np.floor(cell[:, 0]))
    high[:, 0] = np.where(vertical, np.rint(cell[:, 0]), np.ceil(cell[:, 0]))
    low[:, 1] = np.where(vertical, np.floor(cell[:, 1]), np.rint(cell[:, 1]))
    high[:, 1] = np.where(vertical, np.ceil(cell[:, 1]), np.rint(cell[:, 1]))
    low = low * pitch + origin
    high = high * pitch + origin
    fa = solid.evaluate(low[:, 0], low[:, 1], np.asarray(z))
    fb = solid.evaluate(high[:, 0], high[:, 1], np.asarray(z))
    crossing = (fa < 0) != (fb < 0)
    # Exact sampled zeros stay on their original vertices.
    active = crossing & (fa != 0) & (fb != 0)
    for _ in range(8):
        if not active.any():
            break
        middle = (low + high) / 2
        fm = solid.evaluate(middle[:, 0], middle[:, 1], np.asarray(z))
        same = (fm < 0) == (fa < 0)
        move_low = active & same
        move_high = active & ~same
        low[move_low] = middle[move_low]
        fa[move_low] = fm[move_low]
        high[move_high] = middle[move_high]
    result = xy.copy()
    result[active] = ((low + high) / 2)[active]
    return result


@dataclass(frozen=True)
class ImplicitContourSource:
    """Uniform marching squares with root refinement and explicit clipping.

    tolerance limits XY grid pitch. Features smaller than that pitch may be
    missed; field topology and Hausdorff error are not certified. The point cap
    prevents accidentally allocating an enormous plane. No 3D grid is created.
    """

    solid: BoundedSolid
    max_grid_points: int = 4_000_000
    tile_rows: int = 128

    def __post_init__(self) -> None:
        if (
            not isinstance(self.max_grid_points, int)
            or isinstance(self.max_grid_points, bool)
            or not isinstance(self.tile_rows, int)
            or isinstance(self.tile_rows, bool)
            or self.max_grid_points < 4
            or self.tile_rows < 1
        ):
            raise ValueError("Require at least four grid points and positive tile rows")

    def contours_at(self, z: float, *, tolerance: float) -> ContourLayer:
        """Extract closed polygon regions, reporting the sampled-grid limitation."""
        from shapely.geometry import GeometryCollection, Polygon, box
        from skimage.measure import find_contours

        if not math.isfinite(z) or not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("Require finite Z and positive finite contour tolerance")
        position = LayerSpec.plane(z)
        bounds = self.solid.build_bounds
        diagnostics = (
            "Uniform sampling can miss sub-grid features; topology is not certified.",
        )
        if z <= bounds[4] or z >= bounds[5]:
            return ContourLayer(position, (), diagnostics=diagnostics)
        nx = max(2, math.ceil((bounds[1] - bounds[0]) / tolerance) + 1)
        ny = max(2, math.ceil((bounds[3] - bounds[2]) / tolerance) + 1)
        if (nx + 2) * (ny + 2) > self.max_grid_points:
            raise ValueError(
                "Contour grid exceeds max_grid_points; increase tolerance or the cap"
            )
        xs = np.linspace(bounds[0], bounds[1], nx)
        ys = np.linspace(bounds[2], bounds[3], ny)
        pitch = (float(xs[1] - xs[0]), float(ys[1] - ys[0]))
        values = np.empty((ny, nx))
        for start in range(0, ny, self.tile_rows):
            stop = min(start + self.tile_rows, ny)
            values[start:stop] = self.solid.evaluate(
                xs[None, :], ys[start:stop, None], np.asarray(z)
            )
        # A strictly outside halo closes contours for fields touching the frame.
        padded = np.pad(
            values, 1, constant_values=max(1.0, float(np.max(np.abs(values))))
        )
        origin = (bounds[0] - pitch[0], bounds[2] - pitch[1])
        geometry = GeometryCollection()
        for contour in find_contours(padded, 0, fully_connected="low"):
            xy = np.column_stack(
                (
                    origin[0] + contour[:, 1] * pitch[0],
                    origin[1] + contour[:, 0] * pitch[1],
                )
            )
            if len(xy) < 4:
                continue
            xy = _refine_roots(self.solid, xy, pitch=pitch, z=z, origin=origin)
            polygon = Polygon(xy)
            if polygon.area == 0:
                continue
            if not polygon.is_valid:
                raise ValueError(
                    "Contour extraction produced invalid topology; refine the sampling"
                )
            # XOR recovers holes and nested islands without depending on ring order.
            geometry = geometry.symmetric_difference(polygon)
        geometry = geometry.intersection(box(*[bounds[i] for i in (0, 2, 1, 3)]))
        if geometry.is_empty:
            polygons = []
        elif geometry.geom_type == "Polygon":
            polygons = [geometry]
        elif geometry.geom_type == "MultiPolygon":
            polygons = list(geometry.geoms)
        else:
            raise ValueError("Contour extraction produced non-polygon geometry")
        regions = tuple(
            ContourRegion(
                np.asarray(p.exterior.coords),
                tuple(np.asarray(h.coords) for h in p.interiors),
            )
            for p in sorted(polygons, key=lambda polygon: polygon.bounds)
        )
        return ContourLayer(position, regions, pitch, diagnostics)
