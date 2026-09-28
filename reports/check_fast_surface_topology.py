"""Check closure and periodic caps for the selected fast surface settings."""

import json

import graded_comparison as gc
import meshers
import numpy as np
import pyvista as pv
from matched_surface_benchmark import CASES, HERE


def main():
    records = []
    for mode, filename in (
        ("meshers_fast", "fast_surface_data.json"),
        ("meshers_linear", "linear_surface_data.json"),
    ):
        data = json.loads((HERE / filename).read_text(encoding="utf-8"))
        for case, geometry, grading, repeats, thickness in CASES:
            gc.GEOMETRY = geometry
            gc.GRADE_AXES = grading
            gc.UNIFORM_THICKNESS = thickness
            periodic = tuple(axis not in grading for axis in "xyz")
            resolution = data["selected_resolutions"][case]
            surface = meshers.generate_surface(
                gc.normalized_field,
                bounds=(-repeats / 2, repeats / 2) * 3,
                cells=resolution * repeats - 1,
                band=(-1, 1),
                periodic=periodic,
                smoothing_iterations=0,
                improvement_rounds=0,
                polish_passes=0,
                refine_edges=mode != "meshers_linear",
            )
            triangles = surface.triangles
            mesh = pv.PolyData(
                surface.points,
                np.column_stack((np.full(len(triangles), 3), triangles)),
            )
            mismatches = gc.periodic_mismatch(surface.points, triangles, periodic)
            record = {
                "case": case,
                "mode": mode,
                "triangles": len(triangles),
                "open_edges": mesh.n_open_edges,
                **mismatches,
            }
            records.append(record)
            print(json.dumps(record), flush=True)
    (HERE / "fast_surface_topology.json").write_text(
        json.dumps({"runs": records}, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
