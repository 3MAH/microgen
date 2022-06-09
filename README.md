# microgen
Microstructure generation

[![build-and-test workflow](https://github.com/3MAH/microgen/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/3MAH/microgen)
[![Anaconda-Server Badge](https://anaconda.org/set3mah/microgen/badges/installer/conda.svg)](https://conda.anaconda.org/set3mah)
[![PyPI version](https://badge.fury.io/py/microgen.svg)](https://pypi.org/project/microgen/1.0/)
[![3MAH](https://img.shields.io/badge/website-3MAH-blue)](https://3mah.github.io/)


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

With pip:
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
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/shapes.png" width="250" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/basic_shapes.html#platon-polyhedra"> 
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/platon.png" width="250" height="250">
</a>

### Repeated cells

<a href="https://microgen.readthedocs.io/en/latest/examples/lattices.html#octet-truss"> 
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/octettruss.png" width="250" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/lattices.html#honeycomb"> 
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/honeycomb.png" width="250" height="250">
</a>

### Triply Periodic Minimal Surfaces (TPMS)
<a href="https://microgen.readthedocs.io/en/latest/examples/tpms.html#gyroid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/gyroid.png" width="250" height="250"></a>
<a href="https://microgen.readthedocs.io/en/latest/examples/tpms.html#tpms-available">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/tpms.png" width="250" height="250"></a>
<a href="https://microgen.readthedocs.io/en/latest/examples/tpms.html#spherical-gyroid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/tpms_sphere.png" width="250" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/tpms.html#shell">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/tpms_shell.png" width="250" height="250">
</a>

### 3D operations
<a href="https://microgen.readthedocs.io/en/latest/examples/3d_operations.html#repeating-unit-geometry">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/repeated_geometry.png" width="250" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/3d_operations.html#raster-ellipsoid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/raster.png" width="250" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/3d_operations.html#voronoi">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/Voronoi.png" width="250" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/3d_operations.html#voronoi-gyroid">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/voronoi_gyroid.png" width="250" height="250">
</a>

### Mesh
<a href="https://microgen.readthedocs.io/en/latest/examples/mesh.html#id1">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/Mesh.png" width="250" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/mesh.html#periodic-mesh">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/meshPeriodic.png" width="250" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/mesh.html#mmg">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/mmg.png" width="250" height="250">
</a>
<a href="https://microgen.readthedocs.io/en/latest/examples/mesh.html#mmg-voronoi">
    <img src="https://raw.githubusercontent.com/3MAH/microgen/main/docs/_static/mmg-voro.png" width="250" height="250">
</a>