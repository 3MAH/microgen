"""Compare microgen and meshers surfaces at similar triangle counts.

Run with the experimental meshers package on PYTHONPATH and pass its examples
directory with --meshers-examples. Every measurement uses a fresh process.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CASES = (
    ("gyroid_unit", "gyroid", "none", 1, 0.6),
    ("split_p_unit", "split_p", "none", 1, 0.5),
    ("graded_gyroid_2", "gyroid", "xyz", 2, 0.5),
    ("graded_gyroid_3", "gyroid", "xyz", 3, 0.5),
    ("graded_split_p_2", "split_p", "xyz", 2, 0.5),
    ("lateral_gyroid_2", "gyroid", "x", 2, 0.5),
)


def metrics(points, triangles):
    vertices = points[triangles]
    edge_vectors = np.stack(
        (
            vertices[:, 1] - vertices[:, 0],
            vertices[:, 2] - vertices[:, 1],
            vertices[:, 0] - vertices[:, 2],
        ),
        axis=1,
    )
    lengths = np.linalg.norm(edge_vectors, axis=2)
    area = 0.5 * np.linalg.norm(
        np.cross(edge_vectors[:, 0], -edge_vectors[:, 2]), axis=1
    )
    cosines = np.clip(
        (
            lengths[:, :, None] ** 2
            + np.roll(lengths, 1, axis=1)[:, :, None] ** 2
            - np.roll(lengths, -1, axis=1)[:, :, None] ** 2
        )
        / (2 * lengths[:, :, None] * np.roll(lengths, 1, axis=1)[:, :, None]),
        -1,
        1,
    )[:, :, 0]
    minimum_angle = np.min(np.degrees(np.arccos(cosines)), axis=1)
    equivalent_diameter = np.sqrt(4 * area / np.pi)
    return {
        "triangle_count": len(triangles),
        "minimum_angle_degrees": float(np.min(minimum_angle)),
        "angle_p01_degrees": float(np.quantile(minimum_angle, 0.01)),
        "angle_median_degrees": float(np.median(minimum_angle)),
        "fraction_min_angle_below_10_degrees": float(np.mean(minimum_angle < 10)),
        "area_p01": float(np.quantile(area, 0.01)),
        "area_median": float(np.median(area)),
        "area_p99": float(np.quantile(area, 0.99)),
        "area_cv": float(np.std(area) / np.mean(area)),
        "equivalent_diameter_p01": float(np.quantile(equivalent_diameter, 0.01)),
        "equivalent_diameter_median": float(np.median(equivalent_diameter)),
        "equivalent_diameter_p99": float(np.quantile(equivalent_diameter, 0.99)),
    }


def child(case, mode, resolution, validate=False):
    import graded_comparison as gc

    if mode == "vtk":
        from microgen import Tpms
        from microgen.shape.surface_functions import gyroid, split_p
    else:
        import meshers

    name, geometry, grading, repeats, thickness = next(
        row for row in CASES if row[0] == case
    )
    gc.GEOMETRY = geometry
    gc.GRADE_AXES = grading
    gc.UNIFORM_THICKNESS = thickness
    periodic = tuple(axis not in grading for axis in "xyz")
    bounds = (-repeats / 2, repeats / 2) * 3
    start = time.perf_counter()
    if mode == "vtk":
        shape = Tpms(
            gyroid if geometry == "gyroid" else split_p,
            offset=thickness if grading == "none" else gc.thickness,
            repeat_cell=repeats,
            resolution=resolution,
        )
        surface = shape.generate_surface_mesh()
        seconds = time.perf_counter() - start
        points = surface.points
        triangles = surface.faces.reshape(-1, 4)[:, 1:]
        open_edges = surface.n_open_edges
        actual_cells = resolution * repeats - 1
    else:
        surface = meshers.generate_surface(
            gc.normalized_field,
            bounds=bounds,
            cells=resolution * repeats - 1,
            band=(-1, 1),
            periodic=periodic,
            polish_passes=0 if mode != "meshers" else None,
            improvement_rounds=0 if mode != "meshers" else None,
            smoothing_iterations=0 if mode != "meshers" else None,
            refine_edges=mode != "meshers_linear",
        )
        seconds = time.perf_counter() - start
        points = surface.points
        triangles = surface.triangles
        open_edges = None
        actual_cells = surface.diagnostics["background_cells"]
    validation = {}
    if validate:
        import pyvista as pv

        poly = pv.PolyData(
            points, np.column_stack((np.full(len(triangles), 3), triangles))
        )
        validation = {
            "open_edges": poly.n_open_edges,
            **gc.periodic_mismatch(points, triangles, periodic),
        }
    interior = ~np.any(
        np.isclose(np.abs(points), repeats / 2, atol=1e-8, rtol=0), axis=1
    )
    residual = np.abs(np.abs(gc.normalized_field(*points[interior].T)) - 1)
    print(
        json.dumps(
            {
                "case": name,
                "mode": mode,
                "resolution": resolution,
                "actual_cells": actual_cells,
                "generation_seconds": seconds,
                "open_edges": open_edges,
                "implicit_vertex_residual_p95": float(np.quantile(residual, 0.95)),
                **metrics(points, triangles),
                **validation,
            }
        ),
        flush=True,
    )


def launch(case, mode, resolution, examples, env):
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--child",
        case,
        mode,
        str(resolution),
    ]
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
        timeout=180,
        env=env,
    )
    if completed.returncode:
        return {
            "case": case,
            "mode": mode,
            "resolution": resolution,
            "status": "error",
            "reason": completed.stderr[-1200:],
        }
    return {"status": "ok", **json.loads(completed.stdout.splitlines()[-1])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meshers-examples", type=Path)
    parser.add_argument("--child", nargs=3)
    parser.add_argument(
        "--output", type=Path, default=HERE / "matched_surface_data.json"
    )
    args = parser.parse_args()
    if args.child:
        child(args.child[0], args.child[1], int(args.child[2]))
        return
    if args.meshers_examples is None:
        parser.error("--meshers-examples is required")
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        (str(args.meshers_examples.resolve()), env.get("PYTHONPATH", ""))
    )
    data = {
        "measured_at_utc": datetime.now(timezone.utc).isoformat(),
        "target": "VTK at 16 points per cell; tune meshers to closest triangle count",
        "runs": [],
        "selected_resolutions": {},
    }
    for case, *_ in CASES:
        baseline = launch(case, "vtk", 16, args.meshers_examples, env)
        data["runs"].append({"phase": "scan", **baseline})
        print(json.dumps(data["runs"][-1]), flush=True)
        for resolution in (10, 11, 12, 13, 14, 15):
            value = launch(case, "meshers", resolution, args.meshers_examples, env)
            data["runs"].append({"phase": "scan", **value})
            print(json.dumps(data["runs"][-1]), flush=True)
        candidates = [
            run
            for run in data["runs"]
            if run["case"] == case
            and run["mode"] == "meshers"
            and run["status"] == "ok"
        ]
        if not candidates:
            continue
        best = min(
            candidates,
            key=lambda run: abs(
                np.log(run["triangle_count"] / baseline["triangle_count"])
            ),
        )
        data["selected_resolutions"][case] = best["resolution"]
        for trial in range(3):
            for mode in ("vtk", "meshers") if trial % 2 == 0 else ("meshers", "vtk"):
                resolution = 16 if mode == "vtk" else best["resolution"]
                value = launch(case, mode, resolution, args.meshers_examples, env)
                data["runs"].append({"phase": "timed", "trial": trial + 1, **value})
                print(json.dumps(data["runs"][-1]), flush=True)
        args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
