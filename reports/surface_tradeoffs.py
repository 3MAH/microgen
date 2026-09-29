"""Surface tradeoff benchmark. Run with experimental meshers and examples on PYTHONPATH."""

import json
import os
import signal
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from matched_surface_benchmark import CASES, metrics

HERE = Path(__file__).resolve().parent
SELECTED = {
    "gyroid_unit": 12,
    "split_p_unit": 13,
    "graded_gyroid_2": 12,
    "graded_split_p_2": 12,
}
MODES = ("vtk", "mmgs", "fast", "accurate", "quality")


def periodic_metrics(points, triangles, periodic, bounds):
    from scipy.spatial import cKDTree

    nodes, faces, distances, cap_counts = [], [], [], []
    for axis, enabled in enumerate(periodic):
        if not enabled:
            continue
        transverse = [i for i in range(3) if i != axis]
        caps = []
        for plane in bounds[2 * axis : 2 * axis + 2]:
            mask = np.abs(points[:, axis] - plane) < 1e-8
            vertices = points[mask][:, transverse]
            node_keys = set(map(tuple, np.rint(vertices * 1e8).astype(np.int64)))
            cap = triangles[np.all(mask[triangles], axis=1)]
            face_keys = {
                tuple(
                    sorted(
                        map(
                            tuple,
                            np.rint(points[t][:, transverse] * 1e8).astype(np.int64),
                        )
                    )
                )
                for t in cap
            }
            caps.append((vertices, node_keys, face_keys, len(cap)))
        low, high = caps
        nodes.append(len(low[1] ^ high[1]))
        faces.append(len(low[2] ^ high[2]))
        cap_counts.append([low[3], high[3]])
        distances.append(
            float(
                max(
                    cKDTree(low[0]).query(high[0])[0].max(),
                    cKDTree(high[0]).query(low[0])[0].max(),
                )
            )
            if len(low[0]) and len(high[0])
            else None
        )
    b = np.array(bounds)
    violation = np.maximum(np.maximum(b[::2] - points, points - b[1::2]), 0)
    return dict(
        periodic_node_mismatch=nodes,
        periodic_triangle_mismatch=faces,
        periodic_node_distance_max=distances,
        periodic_cap_triangle_counts=cap_counts,
        box_violation_max=float(violation.max()),
    )


def child(case, mode, size, *, return_mesh=False):
    import graded_comparison as gc
    import meshio
    import pyvista as pv
    from microgen import Tpms
    from microgen.external import Mmg
    from microgen.shape.surface_functions import gyroid, split_p

    _, geometry, grading, repeats, thickness = next(r for r in CASES if r[0] == case)
    gc.GEOMETRY, gc.GRADE_AXES, gc.UNIFORM_THICKNESS = geometry, grading, thickness
    periodic = tuple(axis not in grading for axis in "xyz")
    bounds = (-repeats / 2, repeats / 2) * 3
    start = time.perf_counter()
    if mode in ("vtk", "mmgs", "mmgs_count"):
        shape = Tpms(
            gyroid if geometry == "gyroid" else split_p,
            offset=thickness if grading == "none" else gc.thickness,
            repeat_cell=repeats,
            resolution=16,
        )
        surface = shape.generate_surface_mesh()
        points, triangles = surface.points, surface.faces.reshape(-1, 4)[:, 1:]
        if mode.startswith("mmgs"):
            with tempfile.TemporaryDirectory() as tmp:
                source, dest = str(Path(tmp) / "in.mesh"), str(Path(tmp) / "out.mesh")
                meshio.write(source, meshio.Mesh(points, [("triangle", triangles)]))
                Mmg.mmgs(
                    input=source,
                    output=dest,
                    hsiz=size,
                    hausd=0.01 if mode == "mmgs_count" else 0.001,
                    v=-1,
                )
                result = meshio.read(dest)
                points, triangles = result.points, result.cells_dict["triangle"]
        actual_cells = repeats * 16 - 1
    else:
        shape = Tpms(
            gyroid if geometry == "gyroid" else split_p,
            offset=thickness if grading == "none" else gc.thickness,
            repeat_cell=repeats,
            resolution=SELECTED[case],
        )
        surface = shape.generate_meshers_surface(optimization=mode, periodic=periodic)
        points, triangles = surface.points, surface.triangles
        actual_cells = surface.diagnostics["background_cells"]
    elapsed = time.perf_counter() - start
    used, inverse = np.unique(triangles, return_inverse=True)
    unused_points = len(points) - len(used)
    points = points[used]
    triangles = inverse.reshape(-1, 3)
    poly = pv.PolyData(points, np.column_stack((np.full(len(triangles), 3), triangles)))
    vertices = points[triangles]
    caps = np.any(
        np.all(np.isclose(vertices, np.array(bounds)[::2], atol=1e-8, rtol=0), axis=1)
        | np.all(
            np.isclose(vertices, np.array(bounds)[1::2], atol=1e-8, rtol=0), axis=1
        ),
        axis=1,
    )
    wall = vertices[~caps]
    # Sample facet interiors and edge midpoints, not only vertices on exact roots.
    samples = np.concatenate(
        (
            wall.mean(axis=1),
            (wall[:, 0] + wall[:, 1]) / 2,
            (wall[:, 1] + wall[:, 2]) / 2,
            (wall[:, 2] + wall[:, 0]) / 2,
        )
    )
    values = gc.normalized_field(*samples.T)
    step = 1e-6
    gradients = []
    for axis in range(3):
        shift = np.zeros(3)
        shift[axis] = step
        gradients.append(
            (
                gc.normalized_field(*(samples + shift).T)
                - gc.normalized_field(*(samples - shift).T)
            )
            / (2 * step)
        )
    norms = np.linalg.norm(np.array(gradients), axis=0)
    error = np.abs(np.abs(values) - 1) / np.maximum(norms, 1e-12)
    row = dict(
        unused_points=unused_points,
        case=case,
        mode=mode,
        seconds=elapsed,
        hsiz=size if mode.startswith("mmgs") else None,
        actual_cells=actual_cells,
        open_edges=poly.n_open_edges,
        nonmanifold_edges=poly.extract_feature_edges(
            boundary_edges=False,
            non_manifold_edges=True,
            feature_edges=False,
            manifold_edges=False,
        ).n_cells,
        sampled_distance_p95=float(np.quantile(error, 0.95)),
        sampled_distance_max=float(np.max(error)),
        **metrics(points, triangles),
        **periodic_metrics(points, triangles, periodic, bounds),
    )
    if return_mesh:
        return points, triangles, row
    print(json.dumps(row), flush=True)


def launch(case, mode, size=0):
    process = subprocess.Popen(
        [sys.executable, __file__, case, mode, str(size)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=os.environ.copy(),
        start_new_session=os.name != "nt",
    )
    try:
        stdout, stderr = process.communicate(timeout=120)
        if process.returncode:
            return dict(
                case=case, mode=mode, status="error", error=(stderr + stdout)[-2500:]
            )
        return dict(json.loads(stdout.strip().splitlines()[-1]), status="ok")
    except subprocess.TimeoutExpired:
        if os.name == "nt":
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(process.pid)],
                capture_output=True,
                check=False,
            )
        else:
            os.killpg(process.pid, signal.SIGKILL)
        process.communicate()
        return dict(case=case, mode=mode, status="timeout")


def main():
    output = HERE / "surface_tradeoffs.json"
    data = dict(
        method="Three fresh-process runs per path, sequential alternating order. MMGS hsiz tuned separately toward raw VTK triangle count; hausd=0.001. Timings include construction, conversion, MMGS launch and file IO; exclude imports, tuning and validation.",
        runs=[],
        calibration=[],
    )
    for case in SELECTED:
        vtk = launch(case, "vtk")
        if vtk["status"] != "ok":
            raise RuntimeError(vtk)
        target = vtk["triangle_count"]
        # Equilateral edge length at the target average area.
        size = (4 * vtk["area_median"] / np.sqrt(3)) ** 0.5
        candidates = []
        for _ in range(5):
            row = launch(case, "mmgs", size)
            data["calibration"].append(row)
            if row["status"] != "ok":
                break
            candidates.append((abs(row["triangle_count"] / target - 1), size))
            if candidates[-1][0] <= 0.05:
                break
            size *= (row["triangle_count"] / target) ** 0.5
        size = min(candidates)[1] if candidates else size
        for trial in range(3):
            for mode in MODES if trial % 2 == 0 else MODES[::-1]:
                row = launch(case, mode, size)
                row.update(trial=trial + 1, target_triangles=target)
                data["runs"].append(row)
                output.write_text(json.dumps(data, indent=2))
                print(
                    case,
                    mode,
                    row.get("seconds", row.get("error", row["status"])),
                    flush=True,
                )
    for case in SELECTED:
        for mode in MODES:
            rows = [
                r
                for r in data["runs"]
                if r["case"] == case and r["mode"] == mode and r["status"] == "ok"
            ]
            if rows:
                print(
                    case,
                    mode,
                    statistics.median(r["seconds"] for r in rows),
                    rows[0]["triangle_count"],
                    rows[0]["minimum_angle_degrees"],
                )


if __name__ == "__main__":
    if len(sys.argv) == 4:
        child(sys.argv[1], sys.argv[2], float(sys.argv[3]))
    else:
        main()
