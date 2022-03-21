# microgen
Microstructure generation
-------------------------------------------------------------------------------------------------------
To install Microgen with conda: 
```
conda install -c conda-forge -c cadquery -c set3mah microgen
```
You can now use Microgen!

-------------------------------------------------------------------------------------------------------

If you want to modify the sources you would prefer to install microgen from the sources

First, you need to make sure you have all the required dependencies:
```
conda install -c conda-forge -c cadquery -c set3mah python cadquery=2 python-gmsh pygalmesh mmg meshio
```

Then you can install microgen: 
```
python setup.py install
```
