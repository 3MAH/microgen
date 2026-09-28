"""Measure meshers tetrahedra with optimization disabled against raw VTK grids."""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from benchmark_limitations import CASES, MICROGEN_REPO, run_child

SELECTED = {"gyroid_unit", "split_p_unit", "graded_gyroid_2"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--meshers-repo", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("fast_volume_data.json"),
    )
    args = parser.parse_args()
    repo = args.meshers_repo.resolve()
    comparison = repo / "crates/meshers-python/examples/graded_comparison.py"
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        (str(repo / "crates/meshers-python/python"), str(MICROGEN_REPO))
    )
    data = {
        "measured_at_utc": datetime.now(timezone.utc).isoformat(),
        "meshers_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=repo,
            text=True,
        ).strip(),
        "method": "fresh processes, three alternating trials; meshers optimize_passes=0",
        "runs": [],
    }
    for case, geometry, grading, repeats, thickness in CASES:
        if case not in SELECTED:
            continue
        for trial in range(3):
            modes = (
                ("microgen_volume", "meshers_volume")
                if trial % 2 == 0
                else ("meshers_volume", "microgen_volume")
            )
            for mode in modes:
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
                    "--optimize-passes",
                    "0",
                ]
                result = run_child(command, env)
                record = {"case": case, "mode": mode, "trial": trial + 1, **result}
                data["runs"].append(record)
                print(json.dumps(record), flush=True)
        args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
