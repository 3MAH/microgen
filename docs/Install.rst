.. _RST Install:

Install
========

Conda and PyPI packages are available for Microgen

----------------------------------------------------------------------------------

With pip
~~~~~~~~

.. code-block:: bash

   pip install microgen

With conda
~~~~~~~~~~

.. code-block:: bash

   conda install -c conda-forge -c set3mah microgen

----------------------------------------------------------------------------------

To modify the sources, clone this repository and install microgen:

.. code-block:: bash

   git clone https://github.com/3MAH/microgen.git
   cd microgen
   pip install -e .[all]

The `-e` or `--editable` option allows to modify the sources without having to reinstall the package and `[all]` installs the optional development dependencies.
