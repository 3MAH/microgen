"""Replace preliminary direct-native timings with complete microgen API timings."""

import importlib.metadata
import json
import platform
import subprocess
from pathlib import Path
from datetime import datetime, timezone
import meshers
import pyvista
from surface_tradeoffs import HERE, SELECTED, launch

path = HERE / "surface_tradeoffs.json"
data = json.loads(path.read_text())
data["native_only_runs"] = [
    r for r in data["runs"] if r["mode"] in ("fast", "accurate", "quality")
]
data["runs"] = [
    r for r in data["runs"] if r["mode"] not in ("fast", "accurate", "quality")
]
for case in SELECTED:
    target = next(r["target_triangles"] for r in data["runs"] if r["case"] == case)
    for trial in range(3):
        for mode in (
            ("fast", "accurate", "quality")
            if trial % 2 == 0
            else ("quality", "accurate", "fast")
        ):
            row = launch(case, mode)
            row.update(trial=trial + 1, target_triangles=target)
            data["runs"].append(row)
            path.write_text(json.dumps(data, indent=2))
            print(
                case,
                mode,
                row.get("seconds", row.get("error", row["status"])),
                flush=True,
            )
repo = Path(meshers.__file__).resolve().parents[4]
data["metadata"] = dict(
    measured_at=datetime.now(timezone.utc).isoformat(),
    platform=platform.platform(),
    python=platform.python_version(),
    vtk=str(pyvista.vtk_version_info),
    mmgpy=importlib.metadata.version("mmgpy"),
    meshers_source=str(repo),
    meshers_commit=subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repo, text=True
    ).strip(),
    microgen_api="Tpms.generate_meshers_surface; uncommitted experiment changes accompanying this report",
    mmgs="MMGS executable bundled with mmgpy; microgen.external.Mmg.mmgs wrapper; default feature detection",
    caveat="MMGS relaxed and microgen API timings were collected in subsequent sequential blocks, not globally interleaved with all paths.",
)
path.write_text(json.dumps(data, indent=2))
