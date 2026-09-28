"""Collect fresh-process data for the microgen/meshers limitations report.

Run from the microgen experiment checkout with its virtual environment:
    python reports/benchmark_limitations.py --meshers-repo PATH

The meshers checkout must contain the experimental-surfaces native extension.
"""

import argparse
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

MICROGEN_REPO = Path(__file__).resolve().parents[1]
CASES = [
    ("gyroid_unit", "gyroid", "none", 1, 0.6),
    ("split_p_unit", "split_p", "none", 1, 0.5),
    ("graded_gyroid_2", "gyroid", "xyz", 2, 0.5),
    ("graded_gyroid_3", "gyroid", "xyz", 3, 0.5),
    ("graded_split_p_2", "split_p", "xyz", 2, 0.5),
    ("lateral_gyroid_2", "gyroid", "x", 2, 0.5),
]
SURFACE_CASES = {case[0] for case in CASES}
VOLUME_CASES = {"gyroid_unit", "split_p_unit", "graded_gyroid_2", "lateral_gyroid_2"}
CURVED_CASES = ("cylinder_full_wrap", "sphere_full_wrap", "sweep")


def revision(path: Path) -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=path, text=True
    ).strip()


def curved_shape(case: str):
    import numpy as np

    from microgen import CylindricalTpms, SphericalTpms, Sweep
    from microgen.shape.surface_functions import gyroid

    if case == "cylinder_full_wrap":
        return CylindricalTpms(
            radius=1.5,
            surface_function=gyroid,
            offset=0.5,
            resolution=8,
            repeat_cell=(1, 0, 1),
        )
    if case == "sphere_full_wrap":
        return SphericalTpms(
            radius=2.0,
            surface_function=gyroid,
            offset=0.5,
            resolution=8,
            repeat_cell=(1, 0, 0),
        )
    if case == "sweep":
        return Sweep(
            curve_points=np.linspace([0, 0, -1], [0, 0, 1], 16),
            surface_function=gyroid,
            radial_max=1.0,
            offset=0.4,
            resolution=8,
            repeat_cell=(2, 1, 6),
        )
    raise ValueError(case)


def curved_child(case: str, mode: str) -> None:
    start = time.perf_counter()
    shape = curved_shape(case)
    if mode == "microgen_surface":
        mesh = shape.generate_surface_mesh(type_part="sheet")
    elif mode == "microgen_volume":
        mesh = shape.generate_volume_mesh(type_part="sheet")
    else:
        try:
            shape.generate_meshers(type_part="sheet")
        except NotImplementedError as error:
            print(json.dumps({"status": "unsupported", "reason": str(error)}))
            return
        raise AssertionError("Expected this geometry to be unsupported")
    seconds = time.perf_counter() - start
    print(
        json.dumps(
            {
                "status": "ok",
                "runtime_seconds": round(seconds, 6),
                "points": mesh.n_points,
                "elements": mesh.n_cells,
                "open_edges": mesh.n_open_edges if mode == "microgen_surface" else None,
            }
        )
    )


def run_child(command: list[str], env: dict[str, str]) -> dict:
    try:
        completed = subprocess.run(
            command,
            cwd=MICROGEN_REPO,
            env=env,
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "reason": "Exceeded 180 seconds"}
    if completed.returncode:
        return {
            "status": "error",
            "reason": completed.stderr.strip()[-1600:],
        }
    try:
        return {"status": "ok", **json.loads(completed.stdout.splitlines()[-1])}
    except (IndexError, json.JSONDecodeError):
        return {"status": "error", "reason": completed.stdout[-1600:]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--meshers-repo", type=Path)
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument(
        "--group", choices=("all", "surface", "volume", "curved"), default="all"
    )
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("limitations_data.json")
    )
    parser.add_argument("--curved-child", nargs=2, metavar=("CASE", "MODE"))
    args = parser.parse_args()
    if args.curved_child:
        curved_child(*args.curved_child)
        return
    if args.meshers_repo is None:
        parser.error("--meshers-repo is required")
    repo = args.meshers_repo.resolve()
    comparison = repo / "crates/meshers-python/examples/graded_comparison.py"
    extension = repo / "crates/meshers-python/python/meshers/_meshers.pyd"
    if not comparison.is_file() or not extension.is_file():
        parser.error("The meshers checkout needs graded_comparison.py and _meshers.pyd")
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        (str(repo / "crates/meshers-python/python"), str(MICROGEN_REPO))
    )
    result = {
        "measured_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "microgen_commit": revision(MICROGEN_REPO),
        "meshers_commit": revision(repo),
        "trials": args.trials,
        "method": "fresh process per trial; runtime excludes imports and startup",
        "runs": [],
    }
    if args.output.is_file():
        existing = json.loads(args.output.read_text(encoding="utf-8"))
        result["runs"] = existing["runs"]
    if args.group in ("all", "surface", "volume"):
        for case, geometry, grading, repeats, thickness in CASES:
            groups = []
            if case in SURFACE_CASES and args.group in ("all", "surface"):
                groups.append(("surface", ("microgen_surface", "meshers_surface")))
            if case in VOLUME_CASES and args.group in ("all", "volume"):
                groups.append(("volume", ("microgen_volume", "meshers_volume")))
            for kind, modes in groups:
                for trial in range(args.trials):
                    for mode in modes if trial % 2 == 0 else reversed(modes):
                        command = [
                            sys.executable,
                            str(comparison),
                            mode,
                            "--geometry",
                            geometry,
                            "--grade-axes",
                            grading,
                            "--repeats",
                            str(repeats),
                            "--grid-points-per-cell",
                            "16",
                            "--uniform-thickness",
                            str(thickness),
                        ]
                        value = run_child(command, env)
                        run = {
                            "group": kind,
                            "case": case,
                            "mode": mode,
                            "trial": trial + 1,
                            **value,
                        }
                        result["runs"].append(run)
                        print(json.dumps(run), flush=True)
                        args.output.write_text(
                            json.dumps(result, indent=2), encoding="utf-8"
                        )
    if args.group in ("all", "curved"):
        for case in CURVED_CASES:
            for mode in ("microgen_surface", "microgen_volume", "meshers_volume"):
                count = args.trials if mode != "meshers_volume" else 1
                for trial in range(count):
                    command = [
                        sys.executable,
                        str(Path(__file__).resolve()),
                        "--curved-child",
                        case,
                        mode,
                    ]
                    value = run_child(command, env)
                    run = {
                        "group": "curved",
                        "case": case,
                        "mode": mode,
                        "trial": trial + 1,
                        **value,
                    }
                    result["runs"].append(run)
                    print(json.dumps(run), flush=True)
                    args.output.write_text(
                        json.dumps(result, indent=2), encoding="utf-8"
                    )


if __name__ == "__main__":
    main()
