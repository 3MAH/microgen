.. _RST Install:

Install
========

Conda and PyPI packages are available for Microgen

----------------------------------------------------------------------------------

With conda: 
~~~~~~~~~~~

.. code-block:: none

   $ conda install -c conda-forge -c cadquery -c set3mah microgen


With pip:
~~~~~~~~~

.. code-block:: none

   $ pip install microgen

You may need to install dependencies mentioned in the requirements.txt file

.. code-block:: none

   $ pip install -r requirements.txt


----------------------------------------------------------------------------------

To modify the sources, clone this repository and set up the following environment:

Create a conda environment with all the required dependencies

.. code-block:: none

   $ conda env create -f environment.yml -n microgen-dev
   $ conda activate microgen-dev


Then install microgen: 

.. code-block:: none

   $ python setup.py install