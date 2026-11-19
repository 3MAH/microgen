# Changelog

All notable changes to microgen are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Each pull request that changes user-visible behaviour adds a line under
`[Unreleased]`, in the matching `Added` / `Changed` / `Deprecated` / `Removed` /
`Fixed` subsection, with its PR number.

## [Unreleased]

Changes since 2.0.0b1, on the way to 2.0.0.

### Added

- `Phase.from_implicit` and `Phase.from_grid` constructors; `Phase` caches its
  sampled grid ([#144]).
- `Phase.translate`, `Phase.rotate` and `Phase.scale` act on the implicit field,
  bounds, iso-value, period, CAD shape, mesh and cached grid together ([#144]).
- Implicit-first Voronoi and strut lattices: convex `Polyhedron` now carries a
  half-space signed-distance field, and lattices expose `generate_implicit()`
  ([#143]).
- `microgen.cad` sub-package: direct OCP (OCCT) backend with the `CadShape`
  wrapper, primitive builders, `import_step` and the `require_cad` install gate
  ([#142]).
- `Phase.from_shape`, `Phase.from_cad`, `Phase.from_mesh`, `Phase.surface_mesh`,
  `Phase.pieces` (`Piece`) and `Phase.center_of_mass` ([#142]).
- `shape_to_cad` bridge converting any implicit `Shape` to a CAD solid, and
  `mesh_to_periodic_shell` sewing periodic triangle meshes (TPMS, spinodoids)
  into closed OCCT shells ([#141]).
- `Rve.from_min_max`, `Rve.grid()`, `Rve.min_point` / `Rve.max_point` and
  per-axis periodicity flags `pbc` ([#140]).
- `generate_surface_mesh(bounds=..., resolution=...)` samples the implicit field
  on request for `Box`, `Sphere` and `Spinodoid` ([#139]).
- `GradedInfill` and `Sweep` are exported from `microgen.shape` ([#139]).
- `Piece`, `CadShape` and `shape_to_cad` are importable from the top-level
  `microgen` package.
- `examples/ImplicitShapes/gyroid_implicit_phase.py`, a tour of the
  implicit-first workflow.

### Changed

- **Breaking:** `Shape.func` is renamed `Shape.field` ([#144]).
- **Breaking:** affine transforms follow the PyVista convention:
  `translate`, `rotate` and `scale` return a new object by default and modify the
  current one only with `inplace=True` ([#144]).
- `Phase.pieces` labels connected regions with 26-connectivity, so convex
  polyhedra are no longer split at sharp vertices ([#143]).
- **Breaking:** `Phase(shape=...)` is replaced by the `Phase.from_*`
  constructors, and the CAD object is accessed through `Phase.cad` instead of
  `Phase.shape` ([#142]).
- **Breaking:** `Rve` is an immutable (frozen) dataclass ([#140]).
- The shared `try_make_solid` helper of `microgen.shape.periodic_shell` upgrades
  sewn shells to solids for both `Tpms` and `Spinodoid` ([#145]).
- The `honeycomb_*` surface functions are annotated with array inputs and
  outputs, like the other TPMS surface functions ([#147]).

### Removed

- `rotate_pv_euler` ([#139]).
- The Jupyter notebooks of `examples/jupyter_notebooks/` (v1 code) and the
  Binder badge.

## [2.0.0b1] - 2026-05-13

First beta of 2.0, covering all changes since 1.3.2 (including the
`2.0.0dev0` and `2.0.0dev1` pre-releases). See the *Migrating from 1.x*
table below.

### Added

- Implicit (F-rep) modelling: `microgen.shape.implicit_ops` (booleans, smooth
  blends, `batch_smooth_union`, `repeat`, `shell`, `blend`, `from_field`, ...)
  and implicit-field support on `Shape`, with `|`, `&`, `-` and `~` operators
  ([#120]).
- Signed-distance fields on `Box`, `Sphere`, `Cylinder`, `Capsule` and
  `Ellipsoid`, so primitives compose directly with implicit operators ([#133]).
- `Spinodoid` shape built on a Gaussian random field (periodic Fourier-mode
  backend) ([#127]).
- Strut lattices: `AbstractLattice` and 12 presets (`BodyCenteredCubic`,
  `FaceCenteredCubic`, `Cubic`, `Octahedron`, `Cuboctahedron`, `Diamond`,
  `OctetTruss`, `RhombicCuboctahedron`, `RhombicDodecahedron`, `TruncatedCube`,
  `TruncatedCuboctahedron`, `TruncatedOctahedron`) ([#99]).
- `orientation` accepts a `scipy.spatial.transform.Rotation` ([#95]).
- `remesh_keeping_boundaries_for_fem`, with a `periodic` flag ([#103]).
- Python 3.14 support for the core install ([#123]).
- `CITATION.cff` and Zenodo metadata ([#125]).
- Quick start and tutorials in the documentation ([#111]).

### Changed

- **Breaking:** CadQuery is no longer a dependency. The core install
  (mesh + implicit fields) needs only VTK/PyVista; CAD features need the `[cad]`
  extra, which installs `cadquery-ocp-novtk` ([#123]).
- **Breaking:** Python 3.10 or newer is required ([#119], [#123]).
- **Breaking:** `generate_vtk()` is renamed `generate_surface_mesh()` on every
  shape ([#134]).
- Surface meshes generated from implicit fields are watertight: they are capped
  at the bounding box ([#137]).
- The mmg remeshing pipeline retries on failure and snaps coordinate drift that
  broke periodicity checks ([#132]).
- mmg 5.8.0 ([#96]).

### Removed

- **Breaking:** all APIs deprecated during 1.x ([#128]), see the migration
  table below.
- Optional dependencies `ipygany`, `pythreejs` and `ipyvtklink` ([#119]).

### Fixed

- `AbstractLattice` surface meshes ignored `periodic=True` ([#118]).
- Remeshing while keeping boundary elements untouched now works regardless of
  periodicity ([#103]).
- Deprecation warnings from recent PyVista releases ([#126], [#128]).

## Migrating from 1.x

| 1.x | 2.0 |
| --- | --- |
| `pip install microgen` brings CadQuery | `pip install microgen[cad]` for CAD features |
| `shape.generate()` returns a CadQuery shape | `shape.generate_cad()` returns a `microgen.cad.CadShape` |
| `shape.generate_vtk()` / `generateVtk()` | `shape.generate_surface_mesh()` |
| `Phase(shape=cq_shape)` | `Phase.from_cad(cad)`, `Phase.from_shape(shape)`, `Phase.from_mesh(mesh)` |
| `phase.shape` | `phase.cad` |
| `phase.getCenterOfMass()` / `centerOfMass` | `phase.center_of_mass` |
| `phase.getInertiaMatrix()` / `inertiaMatrix` | `phase.inertia_matrix` |
| `Box(dim_x=..., dim_y=..., dim_z=...)` | `Box(dim=(...))` |
| `Ellipsoid(a_x=..., a_y=..., a_z=...)` | `Ellipsoid(radii=(...))` |
| `ExtrudedPolygon(listCorners=...)` | `ExtrudedPolygon(list_corners=...)` |
| `Rve(dim_x=..., ...)`, `rve.x_min`, `rve.dx`, ... | `Rve(dim=...)`, `rve.min_point`, `rve.max_point`, `rve.dim` |
| `schwarzP`, `schoenIWP`, `fischerKochS`, ... | `schwarz_p`, `schoen_iwp`, `fischer_koch_s`, ... |
| `rotateEuler`, `fuseShapes`, `cutPhases`, `repeatShape`, `rasterPhase`, ... | `rotate_euler`, `fuse_shapes`, `cut_phases`, `repeat_shape`, `raster_phase`, ... |
| `meshPeriodic(...)`, `mesh(listPhases=..., mshFileVersion=...)` | `mesh_periodic(...)`, `mesh(list_phases=..., msh_file_version=...)` |
| `periodic(...)` | `periodic_split_and_translate(...)` |
| `remesh_keeping_periodicity_for_fem(...)` | `microgen.remesh.remesh_keeping_boundaries_for_fem(..., periodic=True)` |
| `rotate_pv_euler(...)` | `orientation=Rotation.from_euler(...)` or `rotate(...)` |

Releases before 2.0 are listed on the
[GitHub releases page](https://github.com/3MAH/microgen/releases).

[Unreleased]: https://github.com/3MAH/microgen/compare/v2.0.0b1...HEAD
[2.0.0b1]: https://github.com/3MAH/microgen/compare/v1.3.2...v2.0.0b1
[#95]: https://github.com/3MAH/microgen/pull/95
[#96]: https://github.com/3MAH/microgen/pull/96
[#99]: https://github.com/3MAH/microgen/pull/99
[#103]: https://github.com/3MAH/microgen/pull/103
[#111]: https://github.com/3MAH/microgen/pull/111
[#118]: https://github.com/3MAH/microgen/pull/118
[#119]: https://github.com/3MAH/microgen/pull/119
[#120]: https://github.com/3MAH/microgen/pull/120
[#123]: https://github.com/3MAH/microgen/pull/123
[#125]: https://github.com/3MAH/microgen/pull/125
[#126]: https://github.com/3MAH/microgen/pull/126
[#127]: https://github.com/3MAH/microgen/pull/127
[#128]: https://github.com/3MAH/microgen/pull/128
[#132]: https://github.com/3MAH/microgen/pull/132
[#133]: https://github.com/3MAH/microgen/pull/133
[#134]: https://github.com/3MAH/microgen/pull/134
[#137]: https://github.com/3MAH/microgen/pull/137
[#139]: https://github.com/3MAH/microgen/pull/139
[#140]: https://github.com/3MAH/microgen/pull/140
[#141]: https://github.com/3MAH/microgen/pull/141
[#142]: https://github.com/3MAH/microgen/pull/142
[#143]: https://github.com/3MAH/microgen/pull/143
[#144]: https://github.com/3MAH/microgen/pull/144
[#145]: https://github.com/3MAH/microgen/pull/145
[#147]: https://github.com/3MAH/microgen/pull/147
