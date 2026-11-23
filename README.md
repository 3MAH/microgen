<p align="center">
    <img src="https://github.com/3MAH/microgen/blob/main/docs/_static/microgen.png?raw=true" alt="Microgen logo" width="100%"/>
</p>

**microgen** is a Python library designed to facilitate microstructure generation and meshing. Here are its core features:

- **Triply Periodic Minimal Surfaces (TPMS)**: gyroid, Schwarz, Schoen, Neovius, Fischer-Koch… as sheets or skeletons, graded in density, cell size or cell type, in Cartesian, cylindrical or spherical frames, and as infill of any mesh.
- **Spinodoids**: random, anisotropic bicontinuous microstructures from Gaussian random fields.
- **Lattices**: octet-truss, honeycomb, preset and custom (e.g. auxetic) strut lattices.
- **Virtual composites**: periodic RVEs of spheres, short fibres, ellipsoids and other inclusions.
- **3D Voronoi tessellation**: granular materials and polycrystalline metals, also combined with TPMS.
- **Implicit workflow**: shapes as implicit fields (F-rep), with boolean operations and blending.
- **Meshing**: regular and **periodic** meshing using Gmsh, level-set **remeshing** using Mmg.

Generation of 3D objects is achievable through functions that utilize Open CASCADE (via [cadquery-ocp-novtk](https://github.com/CadQuery/ocp-build-system), the direct OCCT Python binding) or VTK (using PyVista). Neper offers tools for 3D tessellation, while Gmsh handles the generation of both regular and periodic meshes, with Mmg handling remeshing tasks.

The CAD path (`.generate_cad()` on shapes, `Phase`, `fuse_shapes`, periodic split, lattice CAD) is **optional** — install with `pip install 'microgen[cad]'`. The default `pip install microgen` gives you mesh + implicit-field (F-rep) workflows only, which is sufficient for many applications and has a much lighter dependency footprint (no OCCT, no VTK version pin).

<p align="center">
    <a href="https://3mah.github.io/microgen-docs/examples/index.html"><img src="docs/_static/microgen_showcase.png" alt="TPMS, spinodoids and strut lattices generated with microgen" title="TPMS, spinodoids and strut lattices" width="49%"></a>
    <a href="https://3mah.github.io/microgen-docs/examples/Spinodoids/spinodoid.html"><img src="docs/_static/spinodoid_morph.gif" alt="Wireframe spinodoid morphing between two random realisations" title="Spinodoid" width="49%"></a>
</p>

## Hello, gyroid (implicit-first, no CAD needed)

```python
from microgen import Phase, Tpms, surface_functions

gyroid = Tpms(surface_function=surface_functions.gyroid, offset=0.3)
phase  = Phase.from_shape(gyroid, resolution=80)

phase.surface_mesh().save("gyroid.stl")  # marching cubes from the SDF
print(f"COM = {phase.center_of_mass}, pieces = {len(phase.pieces)}")
```

A `Phase` built via `Phase.from_shape` carries the implicit field directly —
marching cubes, connected components (`phase.pieces`), and grid-quadrature
moments all run from the SDF, no OCCT required. CAD (`phase.cad`, STEP export,
gmsh `setPeriodic`) is materialised lazily and only needs the `[cad]` extra.

|        |            |
| -------------- | ---------------- |
| PyPI package   | [![PyPI](https://badge.fury.io/py/microgen.svg)](https://pypi.org/project/microgen/) |
| Conda forge package  | [![Conda](https://anaconda.org/conda-forge/microgen/badges/version.svg)](https://anaconda.org/conda-forge/microgen) |
| Documentation  | [![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://3mah.github.io/microgen-docs/) |
| Status         | [![Status](https://github.com/3MAH/microgen/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/3MAH/microgen) |
| Citation       | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6793573.svg)](https://doi.org/10.5281/zenodo.6793573) |
| License        | [![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |
| Website        | [![Website](https://img.shields.io/badge/website-3MAH-blue)](https://3mah.github.io/) |

## Installation

**Core (mesh + implicit-field / F-rep, no CAD kernel)** — lighter, works with VTK 9.4+:

```bash
pip install microgen
```

**With CAD capabilities** (OCCT via `cadquery-ocp-novtk`, enables `.generate_cad()`, `Phase`, `fuse_shapes`, periodic split, lattice CAD):

```bash
pip install 'microgen[cad]'
```

With conda:

```bash
conda install conda-forge::microgen          # core
conda install conda-forge::microgen ocp      # core + CAD
```

**Python version notes:** Python 3.11 to 3.14 is supported on Linux, macOS (Apple Silicon), and Windows, for both core and `[cad]`, on both PyPI and conda-forge.

What you can do in each mode:

|                                       | Core install | `[cad]` install |
| ------------------------------------- | :----------: | :-------------: |
| `Shape.generate_surface_mesh()` (viz) |      ✅      |        ✅        |
| Implicit fields, TPMS F-rep, booleans |      ✅      |        ✅        |
| `Shape.generate_cad()` (OCCT BREP)        |      ❌      |        ✅        |
| `Phase`, `fuse_shapes`, `cut_*`       |      ❌      |        ✅        |
| Periodic split, lattice CAD export    |      ❌      |        ✅        |
| VTK 9.4+                              |      ✅      |        ✅        |

-------------------------------------------------------------------------------------------------------

To modify the sources, clone this repository and install microgen:

```bash
git clone https://github.com/3MAH/microgen.git
cd microgen
pip install -e .[all]
pre-commit install
```

The `-e` or `--editable` option allows to modify the sources without having to reinstall the package and `[all]` installs the optional development dependencies.


Run tests with pytest:

```bash
pytest tests -n auto
```

## Examples

Click on an image to open the example in the [gallery](https://3mah.github.io/microgen-docs/examples/index.html) of **microgen**'s documentation. The images are rendered by the documentation build, so they always show the current version of each example.

### Triply Periodic Minimal Surfaces (TPMS)

<p>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/tpms_surfaces.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_tpms_surfaces_001.png" alt="Available TPMS" title="Available TPMS" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/gyroid.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_gyroid_001.png" alt="Gyroid" title="Gyroid" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/fischer_koch_s.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_fischer_koch_s_001.png" alt="Fischer-Koch S surface" title="Fischer-Koch S surface" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/tpms_shell.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_tpms_shell_001.png" alt="Gyroid in a shell" title="Gyroid in a shell" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/tpms_sphere.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_tpms_sphere_001.png" alt="Spherical gyroid" title="Spherical gyroid" height="200"></a>
</p>

### Graded TPMS

<p>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/graded_density.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_graded_density_001.png" alt="Density grading" title="Density grading" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/graded_cell_size.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_graded_cell_size_001.png" alt="Cell size grading" title="Cell size grading" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/graded_cell_type.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_graded_cell_type_001.png" alt="Morphing between TPMS types" title="Morphing between TPMS types" height="200"></a>
</p>

### Cylindrical and spherical TPMS

<p>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/cylindrical.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_cylindrical_001.png" alt="Cylindrical gyroid" title="Cylindrical gyroid" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/cylindrical_graded.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_cylindrical_graded_001.png" alt="Radially graded cylindrical gyroid" title="Radially graded cylindrical gyroid" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/spherical.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_spherical_001.png" alt="Spherical gyroid shell" title="Spherical gyroid shell" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/coordinate_frames.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_coordinate_frames_001.png" alt="Spherical, cylindrical and swept frames" title="Spherical, cylindrical and swept frames" height="200"></a>
</p>

### TPMS infill

<p>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/infill_bunny.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_infill_bunny_001.png" alt="Gyroid infill of a mesh" title="Gyroid infill of a mesh" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/infill_graded_distance.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_infill_graded_distance_001.png" alt="Infill graded by the distance to the surface" title="Infill graded by the distance to the surface" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/TPMS/infill_gyroid.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_infill_gyroid_001.png" alt="Gyroid infilled by a gyroid" title="Gyroid infilled by a gyroid" height="200"></a>
</p>

### Spinodoids

<p>
<a href="https://3mah.github.io/microgen-docs/examples/Spinodoids/spinodoid.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_spinodoid_001.png" alt="Spinodoid" title="Spinodoid" height="200"></a>
</p>

### Lattices

<p>
<a href="https://3mah.github.io/microgen-docs/examples/Lattices/octet_truss.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_octet_truss_001.png" alt="Periodic octet-truss mesh" title="Periodic octet-truss mesh" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/Lattices/honeycomb.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_honeycomb_001.png" alt="Honeycomb" title="Honeycomb" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/Lattices/preset_lattices.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_preset_lattices_001.png" alt="Preset strut lattices" title="Preset strut lattices" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/Lattices/custom_auxetic_lattice.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_custom_auxetic_lattice_001.png" alt="Custom auxetic lattice" title="Custom auxetic lattice" height="200"></a>
</p>

### Particle and fibre RVEs

<p>
<a href="https://3mah.github.io/microgen-docs/examples/RVE/spheres_rve.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_spheres_rve_001.png" alt="Periodic RVE of spheres" title="Periodic RVE of spheres" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/RVE/fibres_rve.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_fibres_rve_001.png" alt="Periodic RVE of short fibres" title="Periodic RVE of short fibres" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/Fibers/fibers.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_fibers_001.png" alt="Short fibre composite" title="Short fibre composite" height="200"></a>
</p>

### Implicit workflow

<p>
<a href="https://3mah.github.io/microgen-docs/examples/ImplicitShapes/gyroid_implicit_phase.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_gyroid_implicit_phase_001.png" alt="Gyroid as an implicit Phase" title="Gyroid as an implicit Phase" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/ImplicitShapes/implicit_booleans.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_implicit_booleans_001.png" alt="Boolean operations on implicit fields" title="Boolean operations on implicit fields" height="200"></a>
</p>

### Basic shapes and 3D operations

<p>
<a href="https://3mah.github.io/microgen-docs/examples/BasicShapes/shapes.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_shapes_001.png" alt="Basic shapes" title="Basic shapes" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/BasicShapes/platonic_solids.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_platonic_solids_001.png" alt="Platonic solids" title="Platonic solids" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/3Doperations/repeat_shape.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_repeat_shape_001.png" alt="Repeat a unit cell" title="Repeat a unit cell" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/3Doperations/raster_ellipsoid.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_raster_ellipsoid_001.png" alt="Rasterized ellipsoid" title="Rasterized ellipsoid" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/3Doperations/voronoi.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_voronoi_001.png" alt="Voronoi polycrystal" title="Voronoi polycrystal" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/3Doperations/voronoi_gyroid.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_voronoi_gyroid_thumb.png" alt="Voronoi gyroid" title="Voronoi gyroid" height="200"></a>
</p>

### Meshing and remeshing

<p>
<a href="https://3mah.github.io/microgen-docs/examples/Mesh/periodic_gyroid_remesh.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_periodic_gyroid_remesh_thumb.png" alt="Periodic gyroid mesh from a STEP file" title="Periodic gyroid mesh from a STEP file" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/Mesh/remesh.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_remesh_thumb.png" alt="Periodic remeshing" title="Periodic remeshing" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/Mesh/mmg3d.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_mmg3d_thumb.png" alt="Level-set remeshing with mmg" title="Level-set remeshing with mmg" height="200"></a>
<a href="https://3mah.github.io/microgen-docs/examples/Mesh/mmg_voronoi.html"><img src="https://3mah.github.io/microgen-docs/_images/sphx_glr_mmg_voronoi_thumb.png" alt="Level-set remeshing of a Voronoi mesh" title="Level-set remeshing of a Voronoi mesh" height="200"></a>
</p>

## Citation

If you use **microgen** in academic work, please cite it using the concept DOI, which always resolves to the latest archived version on Zenodo:

> [10.5281/zenodo.6793573](https://doi.org/10.5281/zenodo.6793573)

A BibTeX entry suitable for the latest release:

```bibtex
@software{microgen,
  author    = {Marchais, Kevin and
               Chemisky, Yves and
               Legerstee, Yasmin and
               Guevara Garban, Manuel Ricardo},
  title     = {{microgen: microstructure generation and meshing}},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.6793573},
  url       = {https://doi.org/10.5281/zenodo.6793573}
}
```

To cite a specific release instead of the concept, follow the concept DOI to its Zenodo page and pick the version DOI for the release you used. Machine-readable metadata is also available in [`CITATION.cff`](CITATION.cff).

## Authors and contributors

**microgen** is developed by the [3MAH](https://3mah.github.io/) group and external collaborators.

**Authors**

- Kevin Marchais
- Yves Chemisky
- Yasmin Legerstee
- Manuel Ricardo Guevara Garban

**Contributors**

- Louise Le Barbenchon
- Karim Haffar
- Romain d'Esparbès
- Sudeep Kumar Sahoo
