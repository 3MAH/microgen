"""Probe meshers' core chart abilities beyond microgen's adapter guard.

Run with the experimental meshers Python package first on PYTHONPATH.
These are capability probes, not production meshing paths.
"""

import argparse
import json
import time

import meshers
import numpy as np
import pyvista as pv
from benchmark_limitations import curved_shape


def setup(case, resolution):
    shape = curved_shape(case, resolution=resolution)
    half = shape.cell_size * shape.repeat_cell / 2
    bounds = tuple(v for h in half for v in (-float(h), float(h)))
    k = 2 * np.pi / shape.cell_size
    phase = shape.phase_shift

    def raw(x, y, z):
        return shape.surface_function(
            k[0] * (x + phase[0]),
            k[1] * (y + phase[1]),
            k[2] * (z + phase[2]),
        )

    if case in ("cylinder_full_wrap", "cylinder_sector"):
        radius = shape.cylinder_radius

        def mapping(x, y, z):
            return (
                (radius + x) * np.cos(y / radius),
                (radius + x) * np.sin(y / radius),
                z,
            )

        seam_axis = 1
        collapsed_axis = None
    elif case in ("sphere_full_wrap", "sphere_sector"):
        radius = shape.sphere_radius

        def mapping(x, y, z):
            theta, phi = y / radius + np.pi / 2, z / radius
            return (
                (radius + x) * np.sin(theta) * np.cos(phi),
                (radius + x) * np.sin(theta) * np.sin(phi),
                (radius + x) * np.cos(theta),
            )

        seam_axis = 2
        collapsed_axis = 1 if case == "sphere_full_wrap" else None
    else:
        radial_low = bounds[2]
        angular_low = bounds[4]
        angular_width = bounds[5] - bounds[4]

        def mapping(x, y, z):
            radius = y - radial_low
            angle = (z - angular_low) / angular_width * (2 * np.pi)
            return radius * np.cos(angle), radius * np.sin(angle), x

        seam_axis = 2
        collapsed_axis = 1
    return shape, raw, mapping, bounds, seam_axis, collapsed_axis


def mapped_surface(
    case, raw, mapping, bounds, seam_axis, collapsed_axis, cells, periodic_seam
):
    start = time.perf_counter()
    surface = meshers.generate_surface(
        raw,
        bounds=bounds,
        cells=cells,
        band=(-0.25, 0.25),
        periodic=tuple(axis == seam_axis and periodic_seam for axis in range(3)),
        polish_passes=0,
    )
    generation = time.perf_counter() - start
    x, y, z = surface.points.T
    physical = np.column_stack(mapping(x, y, z))
    poly = pv.PolyData(
        physical,
        np.column_stack((np.full(len(surface.triangles), 3), surface.triangles)),
    )
    labels_to_remove = set()
    if case not in ("cylinder_sector", "sphere_sector"):
        labels_to_remove.update((2 + 2 * seam_axis, 3 + 2 * seam_axis))
    if collapsed_axis is not None:
        labels_to_remove.update((2 + 2 * collapsed_axis, 3 + 2 * collapsed_axis))
    keep = ~np.isin(surface.labels, list(labels_to_remove))
    trimmed = pv.PolyData(
        physical,
        np.column_stack((np.full(np.count_nonzero(keep), 3), surface.triangles[keep])),
    ).clean(tolerance=1e-8, absolute=True)
    return {
        "status": "mesh_returned",
        "seconds": round(generation, 3),
        "triangles": len(surface.triangles),
        "open_edges_after_map_before_seam_cleanup": poly.n_open_edges,
        "open_edges_after_seam_cleanup": trimmed.n_open_edges,
        "trimmed_triangles": trimmed.n_cells,
    }


def mapped_volume(
    case,
    shape,
    raw,
    mapping,
    bounds,
    periodic_seam,
    geometry_tolerance,
    minimum_quality,
):
    start = time.perf_counter()
    if periodic_seam:
        if case != "cylinder_full_wrap":
            raise ValueError("The explicit periodic seam probe is for the cylinder")
        offset = float(shape.offset)
        fields = {
            "lower": lambda x, y, z: raw(x, y, z) - 0.5 * offset,
            "upper": lambda x, y, z: -raw(x, y, z) - 0.5 * offset,
        }
        result = meshers.generate_intersection(
            fields,
            bounds=bounds,
            cells=tuple(int(n) - 1 for n in shape.resolution * shape.repeat_cell),
            coordinate_map=mapping,
            periodic=(False, True, False),
            periodic_transforms={1: np.eye(4)},
            minimum_quality=minimum_quality,
            geometry_tolerance=geometry_tolerance,
        )
    elif case in ("cylinder_full_wrap", "sphere_full_wrap"):
        # Call the real adapter implementation after bypassing its support guard.
        result = shape._mesh_curved_part(
            "sheet",
            minimum_quality=minimum_quality,
            geometry_tolerance=geometry_tolerance,
        )
    elif case in ("cylinder_sector", "sphere_sector"):
        result = shape.generate_meshers(
            minimum_quality=minimum_quality,
            geometry_tolerance=geometry_tolerance,
        )
    else:
        offset = float(shape.offset)
        fields = {
            "lower": lambda x, y, z: raw(x, y, z) - 0.5 * offset,
            "upper": lambda x, y, z: -raw(x, y, z) - 0.5 * offset,
        }
        result = meshers.generate_intersection(
            fields,
            bounds=bounds,
            cells=tuple(int(n) - 1 for n in shape.resolution * shape.repeat_cell),
            coordinate_map=mapping,
            minimum_quality=minimum_quality,
            geometry_tolerance=geometry_tolerance,
        )
    output = {
        "status": "mesh_returned",
        "seconds": round(time.perf_counter() - start, 3),
        "points": len(result.points),
        "tetrahedra": len(result.tetrahedra),
        "minimum_mmg_quality": result.diagnostics.get("minimum_mmg_quality"),
        "sampled_surface_error": result.diagnostics.get("sampled_surface_error"),
    }
    if case == "cylinder_full_wrap":
        rounded = np.rint(result.points * 1e8).astype(np.int64)
        output["coincident_physical_points"] = len(rounded) - len(
            np.unique(rounded, axis=0)
        )
        caps = []
        for tag in (3, 4):
            faces = result.surface[result.boundary_tags == tag]
            caps.append({tuple(sorted(map(tuple, rounded[face]))) for face in faces})
        output["seam_cap_triangles"] = [len(cap) for cap in caps]
        output["seam_caps_match"] = caps[0] == caps[1]
        keep = ~np.isin(result.boundary_tags, [3, 4])
        shell = pv.PolyData(
            result.points,
            np.column_stack(
                (
                    np.full(np.count_nonzero(keep), 3),
                    result.surface[keep],
                )
            ),
        ).clean(tolerance=1e-8, absolute=True)
        output["open_edges_after_seam_weld"] = shell.n_open_edges
        grid = pv.UnstructuredGrid(
            {pv.CellType.TETRA: result.tetrahedra}, result.points
        )
        welded = grid.clean(tolerance=1e-8)
        output["welded_volume_points"] = welded.n_points
        output["welded_volume_tetrahedra"] = welded.n_cells
        output["welded_volume_surface_open_edges"] = welded.extract_surface(
            algorithm=None
        ).n_open_edges
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "case",
        choices=(
            "cylinder_full_wrap",
            "cylinder_sector",
            "sphere_full_wrap",
            "sphere_sector",
            "sweep",
        ),
    )
    parser.add_argument("kind", choices=("surface", "volume"))
    parser.add_argument("--cells", type=int, default=16)
    parser.add_argument("--resolution", type=int, default=8)
    parser.add_argument("--geometry-tolerance", type=float, default=0.05)
    parser.add_argument("--minimum-quality", type=float, default=0)
    parser.add_argument("--periodic-seam", action="store_true")
    args = parser.parse_args()
    shape, raw, mapping, bounds, seam_axis, collapsed_axis = setup(
        args.case, args.resolution
    )
    try:
        if args.kind == "surface":
            output = mapped_surface(
                args.case,
                raw,
                mapping,
                bounds,
                seam_axis,
                collapsed_axis,
                args.cells,
                args.periodic_seam,
            )
        else:
            output = mapped_volume(
                args.case,
                shape,
                raw,
                mapping,
                bounds,
                args.periodic_seam,
                args.geometry_tolerance,
                args.minimum_quality,
            )
    except Exception as error:
        output = {
            "status": "rejected",
            "error_type": type(error).__name__,
            "reason": str(error),
        }
    print(json.dumps({"case": args.case, "kind": args.kind, **output}), flush=True)


if __name__ == "__main__":
    main()
