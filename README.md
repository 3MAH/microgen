<p align="center">
    <img src="https://github.com/3MAH/microgen/blob/main/docs/_static/microgen.png?raw=true" alt="Microgen logo" width="100%"/>
</p>

**Microgen** is a Python library designed to facilitate microstructure generation and meshing. Here are its core features:

- **Repeated cells**: Generation of lattice structures such as octet trusses and honeycombs.
- **Triply Periodic Minimal Surfaces (TPMS)**: TPMS-based lattice generation known for favorable physical (mechanical, thermal, ...) properties like low density and large surface area.
- **Virtual composites microstructures**: Generation of basic reinforcement geometries, including spheres, cylinders, ellipsoids, and more.
- **3D Voronoi tessellation**: Simulation of granular materials and polycrystalline metals.
- **Meshing**: Regular and **periodic** meshing using Gmsh, **remeshing** using Mmg.

Generation of 3D objects is achievable through functions that utilize Open CASCADE (via Cadquery) or VTK (using PyVista). Neper offers tools for 3D tessellation, while Gmsh handles the generation of both regular and periodic meshes, with Mmg handling remeshing tasks.

<p align="center">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/gyroid.gif" alt="Gyroid" width="49%"/>
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/fischerKoch.gif" alt="TPMS" width="49%"/>
</p>

|        |            |
| -------------- | ---------------- |
| PyPI package   | [![PyPI](https://badge.fury.io/py/microgen.svg)](https://pypi.org/project/microgen/) |
| Conda forge package  | [![Conda](https://anaconda.org/conda-forge/microgen/badges/version.svg)](https://anaconda.org/conda-forge/microgen) |
| Documentation  | [![Documentation](https://readthedocs.org/projects/microgen/badge/?version=latest)](https://microgen.readthedocs.io/en/latest/?badge=latest) |
| Status         | [![Status](https://github.com/3MAH/microgen/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/3MAH/microgen) |
| Citation       | [![DOI](https://zenodo.org/badge/380437028.svg)](https://zenodo.org/badge/latestdoi/380437028) |
| License        | [![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |
| Website        | [![Website](https://img.shields.io/badge/website-3MAH-blue)](https://3mah.github.io/) |
| Binder         | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/3MAH/microgen/HEAD?urlpath=lab%2Ftree%2Fexamples%2Fjupyter_notebooks) |

## Installation

With pip:

```bash
pip install microgen
```

With conda:

```bash
conda install conda-forge::microgen
```

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

Click on the image to be redirected to the corresponding example on Microgen's documentation

### Basic shapes

<a href="https://microgen.readthedocs.io/en/latest/examples/basic_shapes.html#basic-shapes">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/shapes.png" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/basic_shapes.html#platon-polyhedra">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/platon.png" height="250">
</a>

### Repeated cells

<a href="https://microgen.readthedocs.io/en/latest/examples/lattices.html#octet-truss">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/octettruss.png" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/lattices.html#honeycomb">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/honeycomb.png" height="250">
</a>

### Triply Periodic Minimal Surfaces (TPMS)

<a href="https://microgen.readthedocs.io/en/latest/examples/tpms.html#gyroid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/gyroid.png"height="250"></a>
<a href="https://microgen.readthedocs.io/en/latest/examples/tpms.html#tpms-available">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/tpms.png" height="250"></a>
<a href="https://microgen.readthedocs.io/en/latest/examples/tpms.html#spherical-gyroid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/tpms_sphere.png" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/tpms.html#shell">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/tpms_shell.png" height="250">
</a>

### 3D operations

<a href="https://microgen.readthedocs.io/en/latest/examples/3d_operations.html#repeating-unit-geometry">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/repeated_geometry.png" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/3d_operations.html#raster-ellipsoid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/rasterEllipsoid.png" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/3d_operations.html#voronoi">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/Voronoi.png" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/3d_operations.html#voronoi-gyroid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/Gyroid-voro.png" height="250">
</a>

### Mesh

<a href="https://microgen.readthedocs.io/en/latest/examples/mesh.html#id1">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/Mesh.png" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/mesh.html#periodic-mesh">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/meshPeriodic.png" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/mesh.html#mmg">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/mmg.png" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/mesh.html#mmg-voronoi">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/mmg-voro.png" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/mesh.html#periodic-remeshing">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/examples/remeshed_mesh.png" height="250">
</a>
