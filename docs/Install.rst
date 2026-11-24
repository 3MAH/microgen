.. _RST Install:

Install
=======

microgen supports Python 3.11 to 3.14 on Linux, macOS (Apple Silicon) and
Windows, with PyPI and conda-forge packages.

Core or CAD install
-------------------

The core install covers implicit fields (F-rep), TPMS and surface meshes. The
``[cad]`` extra adds the OCCT kernel (through ``cadquery-ocp-novtk``), needed
for BREP solids, STEP export, Phase CAD operations and periodic gmsh meshes.
The ``[meshers]`` extra adds the Meshers backend of
``generate_volume_mesh(backend="meshers")``, which meshes tetrahedra directly
from the implicit field (wheels for Linux x86-64, macOS and Windows). ``[all]``
installs every optional feature; on Linux aarch64 it skips Meshers.

.. tab-set::

   .. tab-item:: pip

      .. code-block:: bash

         pip install microgen              # core
         pip install "microgen[cad]"       # core + CAD
         pip install "microgen[meshers]"   # core + Meshers volume meshing
         pip install "microgen[all]"       # core + CAD + Meshers + ray tracing

   .. tab-item:: conda

      .. code-block:: bash

         conda install -c conda-forge microgen        # core
         conda install -c conda-forge microgen ocp    # core + CAD

.. list-table::
   :header-rows: 1
   :widths: 50 25 25

   * -
     - Core install
     - ``[cad]`` install
   * - ``Shape.generate_surface_mesh()``
     - yes
     - yes
   * - Implicit fields, TPMS F-rep, booleans
     - yes
     - yes
   * - ``Shape.generate_cad()`` (OCCT BREP)
     - no
     - yes
   * - ``Phase`` CAD operations, ``fuse_shapes``, ``cut_*``
     - no
     - yes
   * - Periodic split, lattice CAD export
     - no
     - yes

Remeshing with mmg needs the mmg executables:
``conda install -c conda-forge mmg``.

Development install
-------------------

To modify the sources, clone the repository and install microgen in editable
mode with every optional feature and the development and documentation tools:

.. code-block:: bash

   git clone https://github.com/3MAH/microgen.git
   cd microgen
   pip install -e ".[all,dev,docs]"
   pre-commit install

Run the tests with ``pytest tests -n auto``. To build this documentation
locally:

.. code-block:: bash

   pip install -e ".[docs,cad]"
   make -C docs html
