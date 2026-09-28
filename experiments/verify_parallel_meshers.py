"""Compare automatic single-mesh workers and parallel microgen jobs."""

import json
import statistics
import time
from pathlib import Path

import numpy as np
from microgen import Tpms, generate_meshers_parallel
from microgen.shape.surface_functions import gyroid

shapes = [Tpms(gyroid, offset=0.6, resolution=16) for _ in range(4)]
options = {"periodic": (True,) * 3}
reference = shapes[0].generate_meshers(threads=1, **options)
records = []
for name, operation in [
    ("single_one_worker", lambda: [shapes[0].generate_meshers(threads=1, **options)]),
    ("single_automatic", lambda: [shapes[0].generate_meshers(**options)]),
    (
        "four_sequential",
        lambda: [s.generate_meshers(threads=1, **options) for s in shapes],
    ),
    ("four_parallel", lambda: generate_meshers_parallel(shapes, **options)),
]:
    times = []
    for _ in range(3):
        start = time.perf_counter()
        results = operation()
        times.append(time.perf_counter() - start)
        for result in results:
            np.testing.assert_array_equal(result.points, reference.points)
            np.testing.assert_array_equal(result.tetrahedra, reference.tetrahedra)
            for actual, expected in zip(
                result.periodic_pairs, reference.periodic_pairs, strict=True
            ):
                np.testing.assert_array_equal(actual, expected)
            assert result.diagnostics["minimum_mmg_quality"] >= 0.1
            assert result.diagnostics["sampled_surface_error"] <= 0.01
    row = {
        "case": name,
        "seconds": statistics.median(times),
        "times": times,
        "meshes": len(results),
        "diagnostics": results[0].diagnostics,
    }
    print(json.dumps(row), flush=True)
    records.append(row)
Path("reports/microgen_parallel_results.json").write_text(json.dumps(records, indent=2))
