# microgen
Microstructure generation

![build-and-test workflow](https://github.com/3MAH/microgen/actions/workflows/build-and-test.yml/badge.svg)
[![Anaconda-Server Badge](https://anaconda.org/set3mah/microgen/badges/installer/conda.svg)](https://conda.anaconda.org/set3mah)


Documentation
--------------
Provider      | Status
--------      | ------
Read the Docs | [![Documentation Status](https://readthedocs.org/projects/microgen/badge/?version=latest)](https://microgen.readthedocs.io/en/latest/?badge=latest)


Installation
------------

-------------------------------------------------------------------------------------------------------
Install Microgen with conda: 
```
conda install -c conda-forge -c cadquery -c set3mah microgen
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
