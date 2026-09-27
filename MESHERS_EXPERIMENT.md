# Direct periodic TPMS meshing with meshers 0.1.0

This branch now changes TPMS volume generation only. `generate_volume_mesh()`
uses meshers for supported TPMS geometries, with no backend selector. Existing
surface meshes, CAD generation, other shapes, phases, lattices, and MMG remeshing
retain their original implementations.

```python
from microgen import Tpms
from microgen.shape.surface_functions import gyroid

shape = Tpms(gyroid, offset=0.5, resolution=16)
mesh = shape.generate_volume_mesh(
    periodic=(True, True, True),
    minimum_quality=0.1,
    geometry_tolerance=0.01,
)
print(mesh.cell_data['MMGQuality'].min())
mesh.save('gyroid.vtu')
```

Quality uses the MMG tetrahedron measure, with 1 for a regular tetrahedron.
The defaults reject nonempty meshes below 0.1 quality or above 0.01 sampled
geometry error in physical coordinate units. These are acceptance gates, not
adaptive refinement, a certified Hausdorff bound, or a guarantee of FEM accuracy.
A quality failure raises `meshers.MeshingError`; it never triggers MMG or a silent
return of the legacy grid. Callers can request stricter limits.

`periodic` explicitly selects matching axes; it defaults to disabled. The field
and grading must match across each requested pair of faces. The generated mesh
has linear tetrahedra, per-cell `Volume` and `MMGQuality`, diagnostic JSON in
`field_data['meshers_diagnostics']`, and node pairs in `periodic_pairs_x/y/z`.
Mapped meshes also retain rigid periodic transformations. `generate_meshers()`
returns the native meshers result, including boundary tags and diagnostics.

Scalar TPMS offsets keep their historical raw-field units. Constant sheets use
meshers bands; variable sheets use separate upper/lower constraints. Density
fitting measures the generated tetrahedra and leaves the caller unchanged if
meshing fails. Distance-based grading is sampled once and interpolated because
its normalization can depend on the whole grid, not an arbitrary callback batch.

Supported fields compile in the meshers wheel. Imported envelopes and sampled
nodal offsets use callbacks. `compile=False` remains available. Construction
`resolution` retains its per-cell meaning. A meshing-call `resolution=(nx,ny,nz)`
overrides total grid-point counts for the whole domain, including repeats.

## Measured support and remaining gaps

The reproducible script is `experiments/verify_tpms_meshers.py`. The final-default
measurements are in `experiments/tpms_meshers_results.json`; earlier callback-path
measurements are preserved separately. The script records failures as well as
successes, independently checks signed element volumes and quality, and compares
opposite boundary triangle connectivity after applying the periodic node map.
It also verifies physical periodic transformations and boundary closure.

The tested Cartesian sheet cases use offset 0.5 and all three periodic axes.
Fifteen of the sixteen built-in functions passed the 0.1 quality and 0.01 error
limits at a tested resolution. Successful representative cases also cover gyroid
skeletals, thin sheets, density fitting, full density, repeated cells, callable
periodic grading, nodal offsets, infill, and cylindrical/spherical sectors.
A cylindrical sector passed with angular rotation and axial translation pairing.
These examples establish coverage of particular inputs, not every parameter set.

| Feature or case | Evidence and current behavior |
| --- | --- |
| `split_p` at offset 0.5 | Failed minimum quality 0.1 at tested resolutions and optimization counts. The API rejects the mesh. This is a quality limitation for these inputs, not missing field support. |
| Anisotropic repeated gyroid | The tested `(0.5, 1.5, 1)` cell with repeats `(2, 1, 1)` and a phase shift failed the quality gate. Increasing optimization and adjusting grid spacing were also tested; see the recorded attempts. Do not claim all anisotropic configurations fail. |
| Grading incompatible with requested periodic axes | Correctly rejected. The nonperiodic linear-grading case passed with periodicity disabled. Distance grading based on a triangulated envelope also failed periodic matching in the tested case and passed without periodic constraints. |
| Full cylindrical wrap | The direct mapped probe produced acceptable tetrahedra but retained coincident, unjoined seam points. Seam welding or an explicit solver equivalence treatment is unfinished. This is an integration gap, not proof that meshers cannot mesh cylinders. |
| Spherical poles or collapsed radial axes | The direct full-sphere map failed with an inverted background tetrahedron. Regular sectors avoid the singularity and passed. A different chart or singularity treatment is needed. |
| Arbitrary `Sweep` | No validated meshers mapping is implemented here. The existing parametric grid is retained. |
| Graded infill | A tested case failed the 0.01 sampled geometry limit at resolution 16; a finer callback run exceeded the 90-second probe limit. Broader quality/performance support remains unverified. |
| Large background grids | Meshers 0.1.0 accepts 4-128 cells per axis, so this integration requires 5-129 grid points per axis after repeats. The default tetrahedron budget is also finite and configurable. |
| Open zero-thickness TPMS surfaces | Kept as the existing surface API. A solid tetrahedral mesh is a different object. |

For full wraps, poles/collapsed axes, and sweeps, `generate_volume_mesh()` retains
the legacy clipped grid and emits an explicit warning that meshers quality and
periodicity guarantees do not apply. Passing meshers controls to such a chart
raises `NotImplementedError`, as does requesting its native `generate_meshers()`
result. There is no user-selectable backend.

MMG remains available for remeshing existing meshes. This branch neither removes
that dependency nor replaces its boundary-preserving adaptation API.

## Verification and release status

Verified on Windows with Python 3.12.13 and the published meshers 0.1.0 wheel:

- 119 selected regression tests passed with the final compiled defaults.
- Five existing TPMS compatibility checks passed.
- Eleven TPMS integration checks cover periodic triangles, rigid placement,
  positive volume, quality/error acceptance, diagnostic export, failure behavior,
  density fitting, and preservation of legacy surface output.
- Two Spinodoid CAD tests were excluded from the final regression run only after
  reproducing their identical failures on untouched microgen commit `38e1912`.
- Ruff 0.15.12 and `git diff --check` pass.

The full upstream suite, other platforms, FEM solution comparisons, and a
convergence study remain release work. The quality threshold alone does not
establish solution accuracy. Conda packaging has not been updated in this experiment.

```sh
python -m pip install -e '.[dev,cad]' 'cadquery-ocp-novtk<8' 'ruff==0.15.12'
python -m pytest tests/test_meshers_backend.py -q
python experiments/verify_tpms_meshers.py
```

The OCCT pin reproduces the tested environment; it is not a dependency change.
