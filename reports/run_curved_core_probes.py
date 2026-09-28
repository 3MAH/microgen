"""Capture direct curved-chart capability probes in fresh Python processes."""

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = [
    (
        "cylinder_sector",
        "volume",
        "--resolution",
        "12",
        "--geometry-tolerance",
        "0.01",
        "--minimum-quality",
        "0.1",
    ),
    (
        "sphere_sector",
        "volume",
        "--resolution",
        "12",
        "--geometry-tolerance",
        "0.01",
        "--minimum-quality",
        "0.1",
    ),
    ("cylinder_sector", "surface", "--cells", "24"),
    ("sphere_sector", "surface", "--cells", "24"),
    (
        "cylinder_full_wrap",
        "volume",
        "--resolution",
        "12",
        "--geometry-tolerance",
        "0.01",
        "--minimum-quality",
        "0.1",
        "--periodic-seam",
    ),
    ("cylinder_full_wrap", "surface", "--cells", "72", "--periodic-seam"),
    (
        "sphere_full_wrap",
        "volume",
        "--resolution",
        "8",
        "--geometry-tolerance",
        "0.01",
        "--minimum-quality",
        "0.1",
    ),
    ("sphere_full_wrap", "surface", "--cells", "48", "--periodic-seam"),
    (
        "sweep",
        "volume",
        "--resolution",
        "8",
        "--geometry-tolerance",
        "0.01",
        "--minimum-quality",
        "0.1",
    ),
    ("sweep", "surface", "--cells", "48", "--periodic-seam"),
]


def main():
    results = []
    for arguments in CASES:
        completed = subprocess.run(
            [sys.executable, str(HERE / "probe_curved_meshers.py"), *arguments],
            capture_output=True,
            text=True,
            check=False,
            timeout=180,
        )
        if completed.returncode:
            result = {
                "case": arguments[0],
                "kind": arguments[1],
                "status": "process_error",
                "reason": completed.stderr[-1600:],
            }
        else:
            result = json.loads(completed.stdout.splitlines()[-1])
        result["arguments"] = list(arguments[2:])
        results.append(result)
        print(json.dumps(result), flush=True)
    (HERE / "curved_core_probe_data.json").write_text(
        json.dumps(
            {
                "measured_at_utc": datetime.now(timezone.utc).isoformat(),
                "runs": results,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
