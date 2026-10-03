"""Slice a bounded gyroid directly; optionally export UVJ with a real printer profile."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from microgen import surface_functions
from microgen.slicing import (
    BoundedSolid,
    ImplicitContourSource,
    LayerPlan,
    PixelGrid,
    read_uvj_profile,
    slice_rasters,
    write_uvj,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        type=Path,
        help="UVJ exported by UVtools with your printer/material settings",
    )
    parser.add_argument("--output", type=Path, default=Path("gyroid.uvj"))
    parser.add_argument("--cell-size", type=float, default=2.0)
    parser.add_argument("--extent", type=float, default=8.0)
    parser.add_argument("--levelset-halfwidth", type=float, default=0.3)
    args = parser.parse_args()
    if min(args.cell_size, args.extent, args.levelset_halfwidth) <= 0:
        parser.error("Geometry parameters must be positive")
    k = 2 * np.pi / args.cell_size

    def sheet(x, y, z):
        # This is an equation-value offset, not a physical wall thickness.
        return (
            np.abs(surface_functions.gyroid(k * x, k * y, k * z))
            - args.levelset_halfwidth
        )

    solid = BoundedSolid(sheet, (0, args.extent, 0, args.extent, 0, args.extent))
    source = ImplicitContourSource(solid)
    section = source.contours_at(args.extent / 2, tolerance=0.05)
    print(
        f"Midplane: {len(section.regions)} solid regions; no mesh or 3D grid constructed"
    )
    if args.profile is None:
        print("Supply --profile to write raster layers for UVtools")
        return
    profile = read_uvj_profile(args.profile)
    size = profile["Properties"]["Size"]
    grid = PixelGrid(
        size["X"],
        size["Y"],
        (size["Millimeter"]["X"] / size["X"], size["Millimeter"]["Y"] / size["Y"]),
    )
    if args.extent > min(*grid.size_mm):
        parser.error("The geometry must fit within the printer display")
    height = size["LayerHeight"]
    # Use complete printer-height slabs, extending the plan beyond the bounded solid.
    z_max = np.ceil(args.extent / height) * height
    plan = LayerPlan.uniform(z_min=0, z_max=z_max, thickness=height)
    write_uvj(
        slice_rasters(solid, plan, grid=grid, supersampling=2),
        args.output,
        profile=profile,
    )
    print(f"Wrote {len(plan.layers)} layers to {args.output}")


if __name__ == "__main__":
    main()
