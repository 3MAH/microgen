# Meshers 0.1.0 experiment

This branch replaces microgen's solid-meshing implementations with the published
`meshers==0.1.0` wheel. It is an experiment based on microgen commit `38e1912`.
It does not depend on unreleased meshers changes.

## What changed

| Microgen path | Experimental implementation |
| --- | --- |
| Bounded `Shape`, boolean combinations, implicit Voronoi phases | Meshers tetrahedra and their closed boundary |
| Box, sphere, capsule, cylinder, ellipsoid, convex polyhedron surfaces | Meshers applied to the existing implicit field |
| Cartesian TPMS sheets | Two-sided bands, preserving raw-field offset units |
| TPMS skeletals | Negative sublevel fields |
| Callable and nodal TPMS grading | Separate upper/lower constraints; nodal offsets are interpolated |
| Cylindrical and spherical sectors | Named constraints in parameter coordinates with a physical coordinate map |
| Infill and graded infill | TPMS constraints intersected with the object envelope |
| Spinodoid surfaces and volumes | Meshers applied to the existing Fourier field |
| Field-backed Phase surfaces and volumes | Meshers with the phase iso-value subtracted |
| First-order lattice FEM meshing | Periodic implicit strut union, with neighboring struts included |

The PyVista return types remain `PolyData` for surfaces and `UnstructuredGrid`
for volumes. Replaced volume paths return only linear tetrahedra. Surface and
volume calls use the same solid definition. Generic shapes cache native results;
the PyVista views copy their arrays so editing a view does not change the cache.

`generate_meshers()` exposes the native result, diagnostics, periodic pairs, and
VTKHDF export. Treat a generic shape's native cached result as read-only. PyVista
meshes retain periodic pairs in field data. Boundary surfaces also retain face
tags and the native node IDs.

## Deliberately retained

Meshers 0.1.0 does not replace CAD construction, STEP input/output, CAD booleans,
arbitrary existing-mesh remeshing, or higher-order Gmsh elements. ExtrudedPolygon
has no implicit representation in microgen and retains its polygon renderer.
Open TPMS zero-level surfaces retain contouring; a closed solid boundary is not
the same object.

Full angular wraps, collapsed radial axes, spherical poles, and arbitrary Sweep
charts retain the existing parametric grid mesher. These need seam identification,
singularity handling, or a separately validated sweep map. Partial cylindrical
and spherical charts use meshers. Requesting a native meshers result for a retained
chart raises `NotImplementedError`; passing meshers options to its legacy generator
raises `TypeError`. There is no retry with VTK after a meshers failure.

Structured grids remain available for sampling, open contours, offset calibration,
and legacy parametric charts. Geometry definitions, grading, density searches,
CAD-dependent lattice radius calibration, and phase quadrature remain in microgen.

## Behavior to review before a PR

- Primitive surfaces are triangulated solid boundaries. Angular tessellation
  arguments such as `theta_resolution`, `phi_resolution`, and box `level` are
  replaced by grid-point `resolution`. This is an intentional API change.
- Resolution counts grid points, so meshers receives `resolution - 1` cells.
  TPMS and spinodoids retain their per-cell resolution and repeat counts.
- The default sampled geometry tolerance is half the largest background-cell
  spacing. Pass `geometry_tolerance` explicitly for a stricter acceptance gate.
  This is not a guaranteed Hausdorff bound or automatic refinement.
- Callbacks default to `compile=False` because microgen fields can call VTK,
  interpolation, and autograd. `compile=True` is available for supported fields.
  This branch makes no speedup claim, and periodic lattice callbacks are costly.
- Thin features still need a resolution study. Geometry failures and budget limits
  propagate. Primitive and TPMS surface meshes can cost substantially more than
  their former visualization-only renderers.
- Use `periodic=(True, True, True)` explicitly for generic shapes, TPMS and
  spinodoids. Lattice volume meshing defaults to periodic. Matching is validated
  by meshers. Angular periodic maps additionally require rigid transforms through
  the meshers `periodic_transforms` option.
- Imported infill envelopes and curved sectors need broader application-level
  validation. No performance or full upstream-suite certification is claimed.

## Try it

```python
import microgen
from microgen.shape.surface_functions import gyroid

shape = microgen.Tpms(gyroid, offset=0.5, resolution=20)
mesh = shape.generate_volume_mesh(periodic=(True, True, True))
mesh.save("gyroid.vtu")
native = shape.generate_meshers(periodic=(True, True, True))
print(native.diagnostics)
```

```sh
python -m pip install -e '.[dev,cad]' 'cadquery-ocp-novtk<8' 'ruff==0.15.12'
python -m pytest tests/test_meshers_backend.py tests/shapes/test_default_mesh_methods.py tests/test_phase.py tests/test_phase_pieces.py tests/test_spinodoid.py tests/test_implicit_voronoi_lattice.py -q -n 4
```

OCCT 8 bindings raised `Bnd_Box::Limits` conversion errors in existing CAD tests
on this Windows machine. OCCT 7.9.3.1.1 passed those tests. The compatibility pin
above is for reproducing this experiment, not a dependency change to microgen.

## Validation on Windows, Python 3.12

The selected regression run covered 126 tests: 125 passed on the first final run.
The remaining lattice-quadrature check passed after restoring the original bounds
for nonperiodic implicit lattices. All 16 new backend checks passed across the
integration runs, covering element orientation, volume, closure, periodic pairing,
rigid placement, variable offsets, density fitting, curved sectors, infill, empty
solids, and retained angular charts. Ruff 0.15.12 and `git diff --check` passed.

The full upstream test suite was not completed. The larger exploratory TPMS run
was stopped; it is not counted as passing. Existing autograd warnings at angular
coordinate singularities remain during TPMS construction, including for charts
whose final mesh is generated in regular parameter coordinates.
