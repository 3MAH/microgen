"""Area-sampled, bidirectional closest-triangle distances to refined TPMS references."""

import json
import numpy as np
import pyvista as pv
import meshers
import graded_comparison as gc
from surface_tradeoffs import CASES, HERE, child


def poly(points, triangles):
    return pv.PolyData(points, np.column_stack((np.full(len(triangles), 3), triangles)))


def samples(points, triangles, n=10000):
    rng = np.random.default_rng(473)
    p = points[triangles]
    area = np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
    chosen = rng.choice(len(p), n, p=area / area.sum())
    uv = rng.random((n, 2))
    uv[uv.sum(axis=1) > 1] = 1 - uv[uv.sum(axis=1) > 1]
    return (
        p[chosen, 0]
        + uv[:, 0, None] * (p[chosen, 1] - p[chosen, 0])
        + uv[:, 1, None] * (p[chosen, 2] - p[chosen, 0])
    )


def distances(source, target):
    points = source
    return np.abs(
        pv.PolyData(points).compute_implicit_distance(target)["implicit_distance"]
    )


benchmark = json.loads((HERE / "surface_tradeoffs.json").read_text())
rows = []
for case in ("gyroid_unit", "split_p_unit", "graded_gyroid_2", "graded_split_p_2"):
    _, geometry, grading, repeats, thickness = next(r for r in CASES if r[0] == case)
    gc.GEOMETRY, gc.GRADE_AXES, gc.UNIFORM_THICKNESS = geometry, grading, thickness
    periodic = tuple(a not in grading for a in "xyz")
    references = []
    for cells in (63, 95):
        ref = meshers.generate_surface(
            gc.normalized_field,
            bounds=(-repeats / 2, repeats / 2) * 3,
            cells=cells,
            band=(-1, 1),
            periodic=periodic,
            polish_passes=0,
            improvement_rounds=0,
            smoothing_iterations=0,
        )
        references.append(
            (
                poly(ref.points, ref.triangles),
                samples(ref.points, ref.triangles),
                len(ref.triangles),
            )
        )
    consistency = np.concatenate(
        (
            distances(references[0][1], references[1][0]),
            distances(references[1][1], references[0][0]),
        )
    )
    check = dict(
        p95=float(np.quantile(consistency, 0.95)),
        max_sampled=float(consistency.max()),
        triangles=[r[2] for r in references],
    )
    for mode in ("vtk", "fast", "accurate", "quality", "mmgs_count"):
        old = next(
            r for r in benchmark["runs"] if r["case"] == case and r["mode"] == mode
        )
        points, triangles, _ = child(case, mode, old["hsiz"] or 0, return_mesh=True)
        candidate = poly(points, triangles)
        forward = distances(samples(points, triangles), references[1][0])
        reverse = distances(references[1][1], candidate)
        d = np.concatenate((forward, reverse))
        row = dict(
            case=case,
            mode=mode,
            triangles=len(triangles),
            p95=float(np.quantile(d, 0.95)),
            p99=float(np.quantile(d, 0.99)),
            max_sampled=float(d.max()),
            reference_check=check,
        )
        rows.append(row)
        (HERE / "replacement_fidelity.json").write_text(json.dumps(rows, indent=2))
        print(json.dumps(row), flush=True)
