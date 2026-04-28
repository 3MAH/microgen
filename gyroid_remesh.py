"""Generate a gyroid TPMS sheet at resolution 30 and export it to STEP.

Uses the new sewn-shell path: the 6 cell-side faces are kept as ``Geom_Plane``
BREP faces (periodicity-friendly) and the interior TPMS triangles are merged
into a single sewn shell with shared edges via :func:`mesh_to_sewn_shell`.
"""

import time

from microgen import Tpms
from microgen.shape.surface_functions import gyroid

t0 = time.perf_counter()
geometry = Tpms(surface_function=gyroid, density=0.30, resolution=30)
shape = geometry.generate(type_part="sheet")
t1 = time.perf_counter()

shape.export_step("gyroid_sheet_res30.step")
t2 = time.perf_counter()

print(f"generate(): {t1 - t0:.2f}s")
print(f"export_step(): {t2 - t1:.2f}s")
print(f"wrote gyroid_sheet_res30.step")
