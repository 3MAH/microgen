<p align="center">
    <img src="https://github.com/3MAH/microgen/blob/main/docs/_static/microgen.png?raw=true" alt="Microgen logo" width="100%"/>
</p>

**microgen** is a Python library designed to facilitate microstructure generation and meshing. Here are its core features:

- **Repeated cells**: Generation of lattice structures such as octet trusses and honeycombs.
- **Triply Periodic Minimal Surfaces (TPMS)**: TPMS-based lattice generation known for favorable physical (mechanical, thermal, ...) properties like low density and large surface area.
- **Virtual composites microstructures**: Generation of basic reinforcement geometries, including spheres, cylinders, ellipsoids, and more.
- **3D Voronoi tessellation**: Simulation of granular materials and polycrystalline metals.
- **Meshing**: Regular and **periodic** meshing using Gmsh, **remeshing** using Mmg.

Generation of 3D objects is achievable through functions that utilize Open CASCADE (via [cadquery-ocp-novtk](https://github.com/CadQuery/ocp-build-system), the direct OCCT Python binding) or VTK (using PyVista). Neper offers tools for 3D tessellation, while Gmsh handles the generation of both regular and periodic meshes, with Mmg handling remeshing tasks.

The CAD path (`.generate_cad()` on shapes, `Phase`, `fuse_shapes`, periodic split, lattice CAD) is **optional** — install with `pip install 'microgen[cad]'`. The default `pip install microgen` gives you mesh + implicit-field (F-rep) workflows only, which is sufficient for many applications and has a much lighter dependency footprint (no OCCT, no VTK version pin).

<p align="center">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/gyroid.gif" alt="Gyroid" width="49%"/>
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/fischerKoch.gif" alt="TPMS" width="49%"/>
</p>

|        |            |
| -------------- | ---------------- |
| PyPI package   | [![PyPI](https://badge.fury.io/py/microgen.svg)](https://pypi.org/project/microgen/) |
| Conda forge package  | [![Conda](https://anaconda.org/conda-forge/microgen/badges/version.svg)](https://anaconda.org/conda-forge/microgen) |
| Documentation  | [![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://3mah.github.io/microgen-docs/) |
| Status         | [![Status](https://github.com/3MAH/microgen/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/3MAH/microgen) |
| Citation       | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6793573.svg)](https://doi.org/10.5281/zenodo.6793573) |
| License        | [![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |
| Website        | [![Website](https://img.shields.io/badge/website-3MAH-blue)](https://3mah.github.io/) |
| Binder         | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/3MAH/microgen/HEAD?urlpath=lab%2Ftree%2Fexamples%2Fjupyter_notebooks) |

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

**Python version notes:** Python 3.10 to 3.14 is supported on Linux, macOS (Apple Silicon), and Windows, for both core and `[cad]`, on both PyPI and conda-forge.

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

Click on the image to be redirected to the corresponding example on **microgen**'s documentation

### Basic shapes

<a href="https://3mah.github.io/microgen-docs/examples/basic_shapes.html#basic-shapes">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/shapes.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/basic_shapes.html#platon-polyhedra">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/platon.png" height="250">
</a>

### Repeated cells

<a href="https://3mah.github.io/microgen-docs/examples/lattices.html#octet-truss">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/octettruss.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/lattices.html#honeycomb">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/honeycomb.png" height="250">
</a>

### Strut-based Lattices

<a href="https://3mah.github.io/microgen-docs/examples/lattices.html#preset-strut-based-lattices">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/lattices.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/lattices.html#custom-strut-based-lattices">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/auxetic_custom_lattice.png" height="250">
</a>

### Triply Periodic Minimal Surfaces (TPMS)

<a href="https://3mah.github.io/microgen-docs/examples/tpms.html#gyroid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/gyroid.png"height="250"></a>
<a href="https://3mah.github.io/microgen-docs/examples/tpms.html#tpms-available">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/tpms.png" height="250"></a>
<a href="https://3mah.github.io/microgen-docs/examples/tpms.html#spherical-gyroid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/tpms_sphere.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/tpms.html#shell">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/tpms_shell.png" height="250">
</a>

### 3D operations

<a href="https://3mah.github.io/microgen-docs/examples/3d_operations.html#repeating-unit-geometry">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/repeated_geometry.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/3d_operations.html#raster-ellipsoid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/rasterEllipsoid.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/3d_operations.html#voronoi">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/Voronoi.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/3d_operations.html#voronoi-gyroid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/Gyroid-voro.png" height="250">
</a>

### Mesh

<a href="https://3mah.github.io/microgen-docs/examples/mesh.html#id1">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/Mesh.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/mesh.html#periodic-mesh">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/meshPeriodic.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/mesh.html#mmg">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/mmg.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/mesh.html#mmg-voronoi">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/mmg-voro.png" height="250">
</a>
<a href="https://3mah.github.io/microgen-docs/examples/mesh.html#periodic-remeshing">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/remeshed_mesh.png" height="250">
</a>

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
