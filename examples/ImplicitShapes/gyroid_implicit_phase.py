"""
Gyroid as an implicit Phase — no CAD detour.
============================================

End-to-end demonstration of the microgen 2.0 implicit-first pipeline.

A :class:`microgen.Tpms` exposes a scalar SDF directly; wrapping it in a
:class:`microgen.Phase` via :meth:`Phase.from_shape` gives access to
:meth:`phase.surface_mesh`, :attr:`phase.pieces`, :attr:`phase.center_of_mass`,
and the immutable transforms :meth:`phase.translated` / :meth:`phase.scaled`,
none of which require the optional ``[cad]`` extra (no OCCT, no STEP, no gmsh).

Run with:

    python examples/ImplicitShapes/gyroid_implicit_phase.py
"""

from pathlib import Path

import numpy as np

from microgen import Phase, Tpms, surface_functions

OUT = Path(__file__).parent

# %% Build an implicit gyroid TPMS — the field is a callable, no BREP built.
gyroid = Tpms(
    surface_function=surface_functions.gyroid,
    offset=0.3,
    cell_size=1.0,
    repeat_cell=1,
    resolution=20,
)

# %% Wrap the implicit Shape in a Phase.  Field-backed: no CAD required.
phase = Phase.from_shape(gyroid, resolution=80)
print(f"Phase bounds: {phase.bounds}")
print(f"Phase period: {phase.period}  (set by Tpms.cell_size)")
print(f"Phase field-backed: {phase.field is not None}")

# %% Marching-cubes surface mesh straight from the SDF.
surface = phase.surface_mesh()
surface.save(str(OUT / "gyroid_phase_surface.stl"))
print(f"Surface mesh: {surface.n_cells} cells -> gyroid_phase_surface.stl")

# %% Moments computed by grid quadrature — no OCCT needed.
print(f"Center of mass:  {phase.center_of_mass}")
print(f"Inertia tensor diagonal:  {np.diag(phase.inertia_matrix)}")

# %% Connected components via Phase.pieces (scipy.ndimage on the SDF grid).
pieces = phase.pieces
print(f"Pieces: {len(pieces)}  (gyroid 'sheet' part is multiply connected)")
for i, piece in enumerate(pieces[:3]):
    print(f"  piece[{i}]: volume={piece.volume:.4f}, com={piece.com}")

# %% Immutable transforms return new Phase objects.
moved = phase.translated((1.0, 0.0, 0.0))
print(f"\nTranslated phase COM: {moved.center_of_mass}  (shifted by +1 in x)")

scaled = phase.scaled(2.0)
print(f"Scaled phase bounds:  {scaled.bounds}  (cell doubled)")

# %% Compose with another implicit Shape via F-rep boolean ops.
from microgen import Sphere

hollowed = gyroid - Sphere(radius=0.3)  # SDF difference, returns a new Shape
hollow_phase = Phase.from_shape(hollowed, bounds=phase.bounds, resolution=80)
print(f"\nHollowed phase pieces: {len(hollow_phase.pieces)}")
hollow_phase.surface_mesh().save(str(OUT / "gyroid_minus_sphere.stl"))
print("Wrote gyroid_minus_sphere.stl")
