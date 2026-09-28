"""Compare the optimized extractor with microgen VTK, including larger graded boxes."""

import argparse
import json
import os
import platform
import statistics
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matched_surface_benchmark as bench

bench.CASES += (
    ("graded_gyroid_5", "gyroid", "xyz", 5, 0.5),
    ("graded_gyroid_7", "gyroid", "xyz", 7, 0.5),
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--child", nargs=3)
    parser.add_argument(
        "--output", type=Path, default=bench.HERE / "table_surface_data.json"
    )
    args = parser.parse_args()
    if args.child:
        bench.child(args.child[0], args.child[1], int(args.child[2]), validate=True)
        return
    resolutions = json.loads((bench.HERE / "linear_surface_data.json").read_text())[
        "selected_resolutions"
    ]
    # Preserve approximately constant sampling per cell when scaling the domain.
    resolutions.update(graded_gyroid_5=11, graded_gyroid_7=11)
    import meshers
    import pyvista as pv

    meshers_root = Path(meshers.__file__).resolve().parents[4]
    data = {
        "measured_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "pyvista": pv.__version__,
        "vtk": pv.vtk_version_info._asdict(),
        "meshers_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=meshers_root, text=True
        ).strip(),
        "meshers_source_status": subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=meshers_root, text=True
        ).strip(),
        "method": "Three fresh processes per path; alternate order; timings exclude imports and validation. Existing matched-count resolutions retained; larger graded gyroids use 11 meshers versus 16 VTK points per unit cell.",
        "selected_resolutions": resolutions,
        "runs": [],
    }
    for case, *_ in bench.CASES:
        for trial in range(3):
            modes = ("vtk", "meshers_fast", "meshers_linear")
            if trial % 2:
                modes = modes[::-1]
            for mode in modes:
                resolution = 16 if mode == "vtk" else resolutions[case]
                result = subprocess.run(
                    [sys.executable, __file__, "--child", case, mode, str(resolution)],
                    capture_output=True,
                    text=True,
                    timeout=180,
                    env=dict(os.environ),
                    check=False,
                )
                record = {
                    "case": case,
                    "mode": mode,
                    "trial": trial + 1,
                    "status": "error",
                }
                if result.returncode:
                    record["reason"] = result.stderr[-2000:]
                else:
                    record.update(
                        json.loads(result.stdout.splitlines()[-1]), status="ok"
                    )
                data["runs"].append(record)
                args.output.write_text(
                    json.dumps(data, indent=2) + "\n", encoding="utf-8"
                )
        records = [r for r in data["runs"] if r["case"] == case]
        print(
            case,
            {
                mode: statistics.median(
                    r["generation_seconds"]
                    for r in records
                    if r["mode"] == mode and r["status"] == "ok"
                )
                if any(r["mode"] == mode and r["status"] == "ok" for r in records)
                else "failed"
                for mode in modes
            },
            flush=True,
        )


if __name__ == "__main__":
    main()
