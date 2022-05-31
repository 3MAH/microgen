Overview
========

Microstructure generation

----------------------------------------------------------------------------------

Install Microgen with conda: 
.. code-block:: none

   $ conda install -c conda-forge -c cadquery -c set3mah microgen

----------------------------------------------------------------------------------

To modify the sources, clone this repository and set up the following environment:

Create a conda environment with all the required dependencies
.. code-block:: none

   $ conda env create -f environment.yml -n microgen-dev
   $ conda activate microgen-dev


Then install microgen: 
.. code-block:: none

   $ python setup.py install