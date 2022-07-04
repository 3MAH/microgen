
[![Microgen Logo](https://github.com/3MAH/microgen/blob/main/docs/_static/microgen.png?raw=true)](https://github.com/3MAH/microgen)

Microstructure generation

[![build-and-test workflow](https://github.com/3MAH/microgen/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/3MAH/microgen)
[![Anaconda-Server Badge](https://anaconda.org/set3mah/microgen/badges/installer/conda.svg)](https://conda.anaconda.org/set3mah)
[![PyPI version](https://badge.fury.io/py/microgen.svg)](https://pypi.org/project/microgen/1.0/)
[![3MAH](https://img.shields.io/badge/website-3MAH-blue)](https://3mah.github.io/)
[![DOI](https://zenodo.org/badge/380437028.svg)](https://zenodo.org/badge/latestdoi/380437028)


## Documentation

Provider      | Status
--------      | ------
Read the Docs | [![Documentation Status](https://readthedocs.org/projects/microgen/badge/?version=latest)](https://microgen.readthedocs.io/en/latest/?badge=latest)


## Installation

-------------------------------------------------------------------------------------------------------
With conda: 
```
conda install -c conda-forge -c cadquery -c set3mah microgen
```

With pip (Linux only):
```
pip install microgen
```

You may need to install dependencies mentioned in the requirements.txt file
```
pip install -r requirements.txt
```
-------------------------------------------------------------------------------------------------------

To modify the sources, clone this repository and set up the following environment:

Create a conda environment with all the required dependencies
```
conda env create -f environment.yml -n microgen-dev
conda activate microgen-dev
```

Then install microgen: 
```
python setup.py install
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
