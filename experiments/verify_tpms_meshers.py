"""Run isolated TPMS checks against the installed meshers release.

Run from the repository root: python experiments/verify_tpms_meshers.py
Results include failures; this is a capability experiment, not a passing-only test.
"""

from __future__ import annotations

import concurrent.futures
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pyvista as pv

from microgen import CylindricalTpms, GradedInfill, Infill, SphericalTpms, Sweep, Tpms
from microgen.shape.tpms_grading import NormedDistance
from microgen.shape import surface_functions as sf


FUNCTIONS = "gyroid schwarz_p schwarz_d neovius schoen_iwp schoen_frd fischer_koch_s pmy honeycomb lidinoid split_p honeycomb_gyroid honeycomb_schwarz_p honeycomb_schwarz_d honeycomb_schoen_iwp honeycomb_lidinoid".split()
CASES = [f"field:{name}" for name in FUNCTIONS] + [
    "upper",
    "lower",
    "anisotropic_repeated",
    "thin",
    "graded_periodic",
    "graded_nonperiodic",
    "density",
    "full_density",
    "nodal_offset",
    "cylindrical_sector",
    "spherical_sector",
    "infill",
    "full_cylinder",
    "full_sphere",
    "sweep",
    "raw_full_cylinder",
    "raw_full_sphere",
    "anisotropic_balanced",
    "repeat_cells",
    "distance_unconstrained",
    "distance_grading",
    "graded_infill",
    "graded_without_periodicity",
    "cylindrical_periodic_sector",
]


def check(name, resolution, optimize_passes=4):
    kwargs = dict(surface_function=sf.gyroid, offset=0.5, resolution=resolution)
    options = dict(
        periodic=(True,) * 3,
        minimum_quality=0.1,
        geometry_tolerance=0.01,
        optimize_passes=optimize_passes,
    )
    part = "sheet"
    constructor = Tpms
    if name.startswith("field:"):
        kwargs["surface_function"] = getattr(sf, name.split(":")[1])
    elif name in ("upper", "lower"):
        part = name + " skeletal"
    elif name in ("anisotropic_repeated", "anisotropic_balanced"):
        kwargs.update(
            cell_size=(0.5, 1.5, 1), repeat_cell=(2, 1, 1), phase_shift=(0.1, 0.2, 0.3)
        )
    elif name == "repeat_cells":
        kwargs["repeat_cell"] = (2, 2, 2)
    elif name == "thin":
        kwargs["offset"] = 0.1
    elif name == "graded_periodic":
        kwargs["offset"] = lambda x, y, z: 0.6 + 0.1 * np.cos(2 * np.pi * x)
    elif name in ("graded_nonperiodic", "graded_without_periodicity"):
        kwargs["offset"] = lambda x, y, z: 0.6 + 0.2 * x
        if name == "graded_without_periodicity":
            options["periodic"] = (False,) * 3
    elif name in ("distance_grading", "distance_unconstrained"):
        kwargs["offset"] = NormedDistance(
            pv.Sphere(radius=2), boundary_offset=0.8, furthest_offset=0.5
        )
    elif name in ("density", "full_density"):
        kwargs.pop("offset")
        kwargs["density"] = 0.3 if name == "density" else 1.0
    elif "cylinder" in name or "cylindrical" in name:
        constructor = CylindricalTpms
        kwargs.update(
            radius=2, repeat_cell=(1, 0 if "full" in name else 2, 1), offset=1
        )
        options["periodic"] = (False,) * 3
        if name == "cylindrical_periodic_sector":
            options["periodic"] = (False, True, True)
        if name.startswith("raw_"):
            kwargs.update(radius=0.8, cell_size=(0.5, 1.0, 1.0))
    elif "sphere" in name or "spherical" in name:
        constructor = SphericalTpms
        kwargs.update(
            radius=2, repeat_cell=(1, 0, 0) if "full" in name else (1, 2, 2), offset=1
        )
        options["periodic"] = (False,) * 3
        if name.startswith("raw_"):
            kwargs.update(radius=0.8, cell_size=(0.5, 1.0, 1.0))
    elif name in ("infill", "graded_infill"):
        constructor = Infill if name == "infill" else GradedInfill
        kwargs.update(
            obj=pv.Sphere(radius=0.8, theta_resolution=12, phi_resolution=12),
            repeat_cell=1,
        )
        if name == "infill":
            kwargs["offset"] = 0.8
        else:
            kwargs.pop("offset")
            kwargs.update(offset_skin=0.8, offset_core=0.4)
        options["periodic"] = (False,) * 3
    elif name == "sweep":
        constructor = Sweep
        kwargs.update(
            curve_points=np.linspace((0, 0, 0), (0, 0, 2), 12), radial_max=0.5
        )
        options["periodic"] = (False,) * 3
    shape = constructor(**kwargs)
    if name == "anisotropic_balanced":
        options["resolution"] = (
            np.rint(
                (resolution - 1)
                * shape.repeat_cell
                * shape.cell_size
                / min(shape.cell_size)
            ).astype(int)
            + 1
        )
    if name == "distance_unconstrained":
        options["periodic"] = (False,) * 3
    if name == "cylindrical_periodic_sector":
        angle = float(shape.cell_size[1] * shape.repeat_cell[1] / shape.cylinder_radius)
        turn = np.eye(4)
        turn[:2, :2] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
        shift = np.eye(4)
        shift[2, 3] = shape.cell_size[2] * shape.repeat_cell[2]
        options["periodic_transforms"] = {1: turn, 2: shift}
    if name == "nodal_offset":
        shape.offset = 0.6 + 0.1 * np.cos(2 * np.pi * shape.grid.points[:, 0])
    start = time.perf_counter()
    if name.startswith("raw_"):
        result = shape._mesh_curved_part(part, **options)
    else:
        result = shape.generate_meshers(part, **options)
    p = result.points[result.tetrahedra]
    determinants = np.linalg.det(p[:, 1:] - p[:, :1])
    edge_sum = sum(
        np.sum((p[:, i] - p[:, j]) ** 2, axis=1)
        for i in range(4)
        for j in range(i + 1, 4)
    )
    quality = np.sqrt(432 * determinants**2 / edge_sum**3)
    assert len(determinants) and np.all(determinants > 0), (
        "nonpositive or empty tetrahedra"
    )
    assert quality.min() >= 0.1 - 1e-10, "quality gate not satisfied"
    verified_faces = []
    for axis, pairs in enumerate(result.periodic_pairs):
        if not options["periodic"][axis]:
            continue
        lookup = dict(pairs.tolist())
        low = result.surface[result.boundary_tags == 2 * axis + 1]
        high = result.surface[result.boundary_tags == 2 * axis + 2]
        mapped = {tuple(sorted(lookup[int(v)] for v in face)) for face in low}
        expected = {tuple(sorted(face)) for face in high}
        assert mapped == expected, "periodic boundary triangulations differ"
        if len(pairs):
            if result.periodic_transforms is not None:
                transform = result.periodic_transforms[axis]
                mapped_points = (
                    result.points[pairs[:, 0]] @ transform[:3, :3].T + transform[:3, 3]
                )
            else:
                shift = np.eye(3)[axis] * shape.cell_size * shape.repeat_cell
                mapped_points = result.points[pairs[:, 0]] + shape.orientation.apply(
                    shift
                )
            assert np.allclose(
                mapped_points, result.points[pairs[:, 1]], atol=1e-10, rtol=0
            ), "periodic coordinates differ"
        verified_faces.append(len(low))
    surface = pv.PolyData(
        result.points,
        np.column_stack((np.full(len(result.surface), 3), result.surface)),
    )
    assert surface.n_open_edges == 0, "open solid boundary"
    return dict(
        status="passed",
        evaluator=result.diagnostics.get("evaluator"),
        seconds=round(time.perf_counter() - start, 3),
        tetrahedra=len(p),
        minimum_quality=float(quality.min()),
        quality_p01=float(np.quantile(quality, 0.01)),
        volume=float(determinants.sum() / 6),
        sampled_error=result.diagnostics.get("sampled_surface_error"),
        background_cells=result.diagnostics.get("background_cells"),
        periodic_face_counts=verified_faces,
        coincident_points=len(result.points)
        - len(np.unique(np.round(result.points, 10), axis=0)),
    )


def run_case(name):
    attempts = []
    for resolution, passes in ((16, 4), (24, 4), (24, 12), (32, 12)):
        try:
            run = subprocess.Popen(
                [
                    sys.executable,
                    __file__,
                    "--case",
                    name,
                    str(resolution),
                    str(passes),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            try:
                stdout, _ = run.communicate(timeout=90)
            except subprocess.TimeoutExpired:
                if os.name == "nt":
                    subprocess.run(
                        ["taskkill", "/PID", str(run.pid), "/T", "/F"],
                        capture_output=True,
                    )
                else:
                    run.kill()
                run.communicate()
                raise
            result = json.loads(stdout.strip().splitlines()[-1])
        except subprocess.TimeoutExpired:
            result = dict(status="timeout", error="exceeded 90 seconds")
        except Exception as exc:
            result = dict(status="runner_error", error=str(exc))
        attempts.append(dict(resolution=resolution, optimize_passes=passes, **result))
        if (
            result["status"] == "passed"
            or result["status"] == "timeout"
            or name in ("full_cylinder", "full_sphere", "sweep", "graded_nonperiodic")
        ):
            break
    return dict(case=name, attempts=attempts)


if __name__ == "__main__":
    if "--case" in sys.argv:
        try:
            output = check(sys.argv[-3], int(sys.argv[-2]), int(sys.argv[-1]))
        except Exception as exc:
            output = dict(
                status="rejected", exception=type(exc).__name__, error=str(exc)
            )
        print(json.dumps(output))
    else:
        os.environ["PYTHONPATH"] = str(Path(__file__).resolve().parents[1])
        selected = (
            sys.argv[2:] if len(sys.argv) > 1 and sys.argv[1] == "--only" else CASES
        )
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            results = []
            for future in concurrent.futures.as_completed(
                [pool.submit(run_case, n) for n in selected]
            ):
                result = future.result()
                results.append(result)
                print(result["case"], result["attempts"][-1]["status"], flush=True)
        path = Path(__file__).with_name("tpms_meshers_results.json")
        if selected != CASES and path.exists():
            results += [
                r
                for r in json.loads(path.read_text())["cases"]
                if r["case"] not in selected
            ]
        path.write_text(
            json.dumps(
                dict(
                    meshers=importlib.metadata.version("meshers"),
                    python=sys.version,
                    platform=sys.platform,
                    cases=sorted(results, key=lambda r: r["case"]),
                ),
                indent=2,
            ),
            encoding="utf-8",
        )
        print(path)
