# Test the Meshers volume backend

This branch adds `generate_volume_mesh(backend="meshers", ...)` to implicit
`Shape`, Cartesian `Tpms`, and `Spinodoid`. The default remains `backend="vtk"`.
Meshers is optional and imports only when selected. Surface meshing, CAD export,
Gmsh and MMG entry points keep their existing behavior.

## Install the tested candidate

Use Python 3.12 and Rust 1.95.0 to build the candidate used for local validation.
Run these commands from this Microgen checkout in an activated virtual environment:

```bash
rustup toolchain install 1.95.0 --profile minimal
RUSTUP_TOOLCHAIN=1.95.0 CARGO_BUILD_JOBS=2 python -m pip install \
  'meshers @ git+https://github.com/kmarchais/meshers.git@ddc376d76b13fc3ad40f635067ad794a2affff09#subdirectory=crates/meshers-python'
python -m pip install -e '.[dev,meshers]'
```

A native wheel built from that commit also works. The `meshers` extra declares
the intended 0.1 SDK compatibility range. Install the candidate explicitly for
this experiment because the version number alone does not identify its fixes.
The source installation above requires access to the Meshers repository.

## Generate a periodic Gyroid mesh

```python
from microgen import Tpms
from microgen.shape.surface_functions import gyroid

shape = Tpms(gyroid, offset=1.1, resolution=24)
mesh = shape.generate_volume_mesh(
    backend="meshers",
    cells=24,
    order=1,
    periodic=(True, True, True),
    geometry_tolerance=0.002,
    minimum_quality=0.1,
    optimize_passes=1,
    max_tetrahedra=250_000,
)
mesh.save("gyroid.vtu")
```

Meshers accepts its native sizing and admission options through this method.
Explicit `cells` overrides the shape's sampling resolution. Otherwise, the
adapter uses the existing grid's interval counts. Its defaults are one thread
and minimum quality 0.1. It makes one meshing call and propagates failures.
There is no automatic MMG step, fallback, retry or relaxation of the guards.
`band` is reserved for the shape definition and cannot be supplied as an option.

TPMS `sheet`, `upper skeletal`, and `lower skeletal` preserve Microgen's raw
surface-function offset units. Existing density calibration supplies the offset.
Scalar offsets and pointwise callable grading are supported. Arrays and
`OffsetGrading` objects raise `NotImplementedError` because their sampled values
need an explicit interpolation policy. Mapped or clipped TPMS subclasses also
raise `NotImplementedError` until they have a domain adapter.

For graded specimens, request periodic axes only where the complete field is
periodic. A linear thickness gradient along z cannot be periodic along z.

## Preserve a Spinodoid specimen

```python
import numpy as np
from microgen import Spinodoid

shape = Spinodoid(
    k0=4 * np.pi, bandwidth=3, offset=0, resolution=16,
    anisotropy=(1, 1, 2), seed=20261026,
    mode_amplitude_threshold=0.05,
)
mesh = shape.generate_volume_mesh(
    backend="meshers", cells=16, periodic=(True, True, True),
    adapt=True, adapt_tolerance=0.0018, geometry_tolerance=0.002,
    min_length=0.002, max_length=0.1, gradation=0.6,
    optimize_passes=1, max_tetrahedra=250_000,
)
```

`Spinodoid.resolution` defines the Fourier realization. Meshing `cells` is
independent. Both backends use the shape's existing field evaluator. Meshers
compiles that callable, including its current modes, amplitudes, phases and
threshold. The adapter preserves the local `[0, L]` origin and then applies the
shape's rotation and center to the mesh. It does not regenerate the random field
or maintain a separate spinodoid formula.
These settings are a tested specimen, not universal defaults for every field.

## Inspect the output

The result is a PyVista `UnstructuredGrid` containing Tet4 or VTK-ordered Tet10.
`MMGQuality` and `CornerSignedVolume` are cell arrays computed from the four
corners. Corner volumes are not curved Tet10 volumes or an FEA certificate.

Field arrays survive a `.vtu` save and reload:

- `meshers_diagnostics` contains the native diagnostics as JSON.
- `meshers_version` records the installed SDK version.
- `meshers_boundary_faces` and `meshers_boundary_tags` retain native faces and
  tags. Tags 1 through 6 identify the low and high caps of local x, y and z.
- `periodic_pairs_x`, `periodic_pairs_y`, and `periodic_pairs_z` contain point-ID
  pairs for requested axes with boundary nodes.
- `periodic_shift_x`, `periodic_shift_y`, and `periodic_shift_z` contain the
  corresponding translation vectors in world coordinates, including rotation.

For a rotated shape, the pair suffix names the local axis. Pair arrays include
Tet10 midside nodes. These arrays do not create solver-specific periodic boundary
conditions automatically.

## Run the integration checks

```bash
python -m pytest tests/test_meshers_backend.py tests/shapes/test_default_mesh_methods.py -q
```

Without Meshers, the optional native tests skip. Default-backend and missing-SDK
tests still run. With the candidate installed, the tests exercise actual native
meshing, exact Fourier-field reproduction, placement, positive corner volumes,
quality admission, periodic node and face bijections, and metadata persistence.
They do not establish a general speed or accuracy advantage over Microgen+MMG.

Local validation used macOS arm64, Python 3.12.9, Meshers candidate `ddc376d`,
NumPy 2.5.4, SciPy 1.18.1, PyVista 0.49.1 and VTK 9.7.1. VTK used its sequential
backend. Two existing CAD spinodoid tests fail in `meshbridge.py` with
`ShellCreationError`. Both failures reproduce on unchanged upstream commit
`38e19120cd4d41fe110528b3bdb59e3c88627dcb` with identical dependencies. The
optional native backend does not change that CAD conversion path.

The local backend and broader implicit regression selection passed 138 tests:

```bash
python -m pytest tests/test_meshers_backend.py \
  tests/shapes/test_default_mesh_methods.py tests/shapes/test_implicit_ops.py \
  tests/shapes/test_tpms_frep.py tests/shapes/test_tpms_grading.py \
  tests/shapes/test_shape_period.py tests/shapes/test_primitive_implicit.py \
  tests/test_implicit_voronoi_lattice.py -q
```

The existing CI does not install the optional SDK. Native meshing coverage
requires installing the candidate and running the command above.

The shared-evaluator follow-up also ran `tests/test_spinodoid.py` with the two
confirmed upstream CAD failures deselected. All 148 selected tests passed.
NumPy evaluation retained bit-identical values against upstream on the six
fixed resolution and point-count controls, using the existing matrix path.
