"""Add MMGS at relaxed tolerance to the surface tradeoff data."""

import json
import math
from surface_tradeoffs import HERE, SELECTED, launch

path = HERE / "surface_tradeoffs.json"
data = json.loads(path.read_text())
data["runs"] = [r for r in data["runs"] if r["mode"] != "mmgs_count"]
for case in SELECTED:
    vtk = next(r for r in data["runs"] if r["case"] == case and r["mode"] == "vtk")
    target = vtk["triangle_count"]
    size = math.sqrt(4 * vtk["area_median"] / math.sqrt(3))
    candidates = []
    for _ in range(8):
        row = launch(case, "mmgs_count", size)
        data["calibration"].append(row)
        path.write_text(json.dumps(data, indent=2))
        if row["status"] != "ok":
            break
        gap = abs(row["triangle_count"] / target - 1)
        candidates.append((gap, size))
        print("calibrate", case, size, row["triangle_count"], flush=True)
        if gap <= 0.05:
            break
        size *= math.sqrt(row["triangle_count"] / target)
    size = min(candidates)[1] if candidates else size
    for trial in range(3):
        row = launch(case, "mmgs_count", size)
        row.update(trial=trial + 1, target_triangles=target)
        data["runs"].append(row)
        path.write_text(json.dumps(data, indent=2))
        print(case, row.get("seconds", row["status"]), flush=True)
