"""Measure unpolished meshers surfaces against raw VTK at similar counts."""

import argparse
import json
import os
import statistics
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from matched_surface_benchmark import CASES, HERE, launch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meshers-examples", type=Path, required=True)
    parser.add_argument(
        "--mode", choices=("meshers_fast", "meshers_linear"), default="meshers_fast"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output is None:
        args.output = HERE / (
            "linear_surface_data.json"
            if args.mode == "meshers_linear"
            else "fast_surface_data.json"
        )
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        (str(args.meshers_examples.resolve()), env.get("PYTHONPATH", ""))
    )
    reference = json.loads(
        (HERE / "matched_surface_data.json").read_text(encoding="utf-8")
    )
    data = {
        "measured_at_utc": datetime.now(timezone.utc).isoformat(),
        "meshers_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=args.meshers_examples.resolve().parents[2],
            text=True,
        ).strip(),
        "method": f"{args.mode}: polish_passes=0, improvement_rounds=0, smoothing_iterations=0; VTK at 16 points per cell",
        "runs": [],
        "selected_resolutions": {},
    }
    for case, *_ in CASES:
        target = next(
            run["triangle_count"]
            for run in reference["runs"]
            if run["phase"] == "timed" and run["case"] == case and run["mode"] == "vtk"
        )
        candidates = []
        for resolution in range(10, 16):
            result = launch(case, args.mode, resolution, args.meshers_examples, env)
            record = {"phase": "scan", **result}
            data["runs"].append(record)
            if result["status"] == "ok":
                candidates.append(result)
            print(json.dumps(record), flush=True)
        if not candidates:
            continue
        best = min(candidates, key=lambda run: abs(run["triangle_count"] / target - 1))
        data["selected_resolutions"][case] = best["resolution"]
        for trial in range(3):
            order = ("vtk", args.mode) if trial % 2 == 0 else (args.mode, "vtk")
            for mode in order:
                resolution = 16 if mode == "vtk" else best["resolution"]
                result = launch(case, mode, resolution, args.meshers_examples, env)
                record = {"phase": "timed", "trial": trial + 1, **result}
                data["runs"].append(record)
                print(json.dumps(record), flush=True)
        args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
        trial_runs = [
            run
            for run in data["runs"]
            if run["case"] == case and run["phase"] == "timed"
        ]
        vtk = statistics.median(
            run["generation_seconds"] for run in trial_runs if run["mode"] == "vtk"
        )
        fast = statistics.median(
            run["generation_seconds"] for run in trial_runs if run["mode"] == args.mode
        )
        print(
            f"SUMMARY {case}: VTK {vtk:.3f}s, {args.mode} {fast:.3f}s, ratio {fast / vtk:.2f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
