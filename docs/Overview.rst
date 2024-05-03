.. _RST Overview:

Overview
========

Microgen is a simple python library that helps to generate and mesh microstructures.

Here are the main features:
- Entirely written in Python 3 and based on `CadQuery <https://cadquery.readthedocs.io/en/latest/>`_ and `PyVista <https://docs.pyvista.org/index.html>`_ libraries.
- It allows to generate simple reinforcement geometries (spheres, cylinder, ellipsoids, …) to generate virtual composites microstructures.
- Three-dimensional Voronoi tessellation rallons to simulate the response granular materials and polycrystalline metals.
- Regular mesh and periodic mesh are implemented using `Gmsh <https://gmsh.info/>`_, remeshing using `Mmg <https://www.mmgtools.org/>`_ is also implemented.


.. image:: https://anaconda.org/set3mah/microgen/badges/version.svg
    :target: https://anaconda.org/set3MAH/microgen

.. image:: https://badge.fury.io/py/microgen.svg
    :target: https://pypi.org/project/microgen/

.. |3MAH| image:: https://3mah.github.io/assets/images/logo_3mah/3mah_logo_vsmall.png
    :width: 50
    :alt: 3MAH website
    :target: https://3mah.github.io/

.. |GitHub| image:: https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png
    :width: 50
    :alt: GitHub repository
    :target: https://github.com/3MAH/microgen

-------------------------------------------------------------------------

 +----------+-----------------------------------------------------------+
 | |3MAH|   |  `3MAH website <https://3mah.github.io/>`_                |
 +----------+-----------------------------------------------------------+
 | |GitHub| |  `GitHub repository <https://github.com/3MAH/microgen>`_  |
 +----------+-----------------------------------------------------------+

-------------------------------------------------------------------------

.. raw:: html

   <p align="center">
      <img src="_static/gyroid.gif" alt="Gyroid" width="49%"/>
      <img src="_static/fischerKoch.gif" alt="TPMS" width="49%"/>
   </p>

Brief examples
--------------

.. jupyter-execute::
   :hide-code:

   import pyvista
   pyvista.set_jupyter_backend('static')
   pyvista.global_theme.background = 'white'
   pyvista.global_theme.window_size = [600, 400]
   pyvista.global_theme.axes.show = False
   pyvista.global_theme.smooth_shading = True
   pyvista.global_theme.split_sharp_edges = True


.. jupyter-execute::

   import microgen

   geometry = microgen.Tpms(
      surface_function=microgen.surface_functions.gyroid,
      offset=0.3
   )
   shape = geometry.sheet

   shape.plot(color='white')
