# microgen
Microstructure generation
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
