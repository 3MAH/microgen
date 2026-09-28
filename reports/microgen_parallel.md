# Parallel meshing in the microgen experiment

Native meshing now defaults to automatic CPU selection through `threads=None`.
This applies to scalar fields and intersections routed through `microgen._meshers`.
It does not change the retained VTK surface workflow or the browser WASM build.

For independent meshes, `generate_meshers_parallel` runs jobs concurrently and
returns native meshers results in input order. Its default is one native worker
per job, up to the available CPU budget. Compiled fields release the Python GIL
during native meshing. Python callback fields can still serialize on the GIL.

```python
from microgen import Tpms, generate_meshers_parallel
from microgen.shape.surface_functions import gyroid

shapes = [Tpms(gyroid, offset=offset, resolution=16)
          for offset in (0.5, 0.6, 0.7, 0.8)]
meshes = generate_meshers_parallel(shapes, periodic=(True, True, True))

# Limit simultaneous jobs or assign more native workers per job.
meshes = generate_meshers_parallel(shapes, max_workers=2, threads=2)
```

`threads=None` on a batch divides its CPU budget across the scheduled jobs.
An explicit worker count reduces the maximum number of concurrent jobs so their
combined worker count stays within that budget. Separate simultaneous batch
calls have separate budgets. Memory use also grows with concurrent jobs; use
`max_workers` to limit large batches.

Pass distinct shape instances and do not mutate shared field state during a
batch. Density fitting may change a shape's offset. A failed job raises instead
of returning a partial list; running jobs finish before the call returns.

## Measured results

Measured on 2026-09-28 with a budget of 24 logical CPUs, using the native release
build from meshers experimental commit `4cf0495`. These results do not describe
the published 0.1.0 wheel. Each timing is the median of three warm runs through
the microgen API, excluding shape construction and result verification.

The case is a periodic gyroid sheet, offset 0.6, resolution 16, with 15,909
tetrahedra per mesh. Four-job tests generate four identical cases using distinct
shape instances.

| Work | Workers | Time |
|---|---|---:|
| One mesh | 1 native | 0.890 s |
| One mesh | 24 native, automatic | 1.541 s |
| Four meshes sequentially | 1 native per job | 3.604 s |
| Four meshes concurrently | 4 jobs, 1 native each | 0.877 s |

The batch achieved 4.11 times the sequential throughput. Automatic worker
selection made this small single mesh 1.73 times slower. More CPUs are not
automatically faster; `threads=1` remains useful for small individual meshes.
Larger meshes and other geometries need separate measurements before selecting
an adaptive worker policy.

Every run produced exactly equal point coordinates, tetrahedra and periodic
pairings. Minimum MMG quality was 0.14404, with no elements below 0.1; sampled
surface error was 0.00340 against a 0.01 tolerance. All fields compiled with zero
Python callbacks. Parallelism did not trade mesh quality for speed in this test.

Raw timings: [microgen_parallel_results.json](microgen_parallel_results.json).
Reproduce with `python experiments/verify_parallel_meshers.py` from this clone,
using the experimental meshers Python package on `PYTHONPATH`.
