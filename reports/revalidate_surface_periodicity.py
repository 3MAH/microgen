import json
from surface_tradeoffs import HERE, launch

path = HERE / "surface_tradeoffs.json"
data = json.loads(path.read_text())
data.setdefault("periodicity_revalidation", [])
for case in ("gyroid_unit", "split_p_unit"):
    for mode in ("vtk", "mmgs", "mmgs_count", "fast", "accurate", "quality"):
        rows = [r for r in data["runs"] if r["case"] == case and r["mode"] == mode]
        if all("periodic_cap_triangle_counts" in row for row in rows):
            continue
        check = launch(case, mode, rows[0]["hsiz"] or 0)
        assert check["status"] == "ok", check
        assert check["triangle_count"] == rows[0]["triangle_count"]
        for row in rows:
            for key in (
                "periodic_node_mismatch",
                "periodic_triangle_mismatch",
                "unused_points",
                "periodic_node_distance_max",
                "periodic_cap_triangle_counts",
                "box_violation_max",
            ):
                row[key] = check[key]
        data["periodicity_revalidation"].append(check)
        print(
            case,
            mode,
            check["unused_points"],
            check["periodic_node_mismatch"],
            check["periodic_triangle_mismatch"],
            check["periodic_node_distance_max"],
            flush=True,
        )
path.write_text(json.dumps(data, indent=2))
