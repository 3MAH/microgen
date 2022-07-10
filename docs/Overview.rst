.. _RST Overview:

Overview
========

Microstructure generation tool based on `CadQuery <https://cadquery.readthedocs.io/en/latest/>`_ library.


.. image:: https://anaconda.org/set3mah/microgen/badges/installer/conda.svg
    :target: https://conda.anaconda.org/set3mah/
    
.. image:: https://badge.fury.io/py/microgen.svg
    :target: https://pypi.org/project/microgen/1.0/

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


Table of contents
-----------------

Introduction
~~~~~~~~~~~~

   - :ref:`RST Overview`
   - :ref:`RST Install`

Documentation
~~~~~~~~~~~~~

   - :ref:`RST Doc`

Examples
~~~~~~~~

   - :ref:`RST Geometries`
   - :ref:`RST Lattices`
   - :ref:`RST Tpms`
   - :ref:`RST 3Dop`
   - :ref:`RST Mesh`




Brief examples
--------------

.. jupyter-execute::
   :hide-code:

   import pyvista
   pyvista.set_jupyter_backend('pythreejs')
   pyvista.global_theme.background = 'white'
   pyvista.global_theme.window_size = [690, 400]
   pyvista.global_theme.axes.show = False
   pyvista.global_theme.smooth_shading = True
   pyvista.global_theme.antialiasing = True

.. jupyter-execute::
   
   import microgen

   geometry = microgen.Tpms(
      surface_function=microgen.tpms.gyroid,
      type_part="sheet",
      thickness=0.05
   )
   shape = geometry.generateVtk()

   shape.plot(color='white')

.. jupyter-execute::

   import cadquery as cq

   capsule = microgen.Capsule(center=(0, 0, 0.5), height=3, radius=1)
   shape = capsule.generate()

   shell = cq.Workplane().add(shape).shell(0.025).split(keepBottom=True).val()
   half_capsule = cq.Workplane().add(shape).split(keepBottom=True).val()

   gyroid = microgen.Tpms(center=(0., 0., 0),
                surface_function=microgen.tpms.gyroid,
                type_part='sheet',
                thickness=0.05,
                cell_size=1,
                repeat_cell=(5, 3, 1))
   shape_gyroid = gyroid.generate()

   inner_gyroid = shape_gyroid.intersect(half_capsule)
   fuse = inner_gyroid.fuse(shell)

   mesh = pyvista.PolyData(fuse.toVtkPolyData(0.1))
   pl = pyvista.Plotter()
   pl.add_mesh(mesh, color='white')
   pl.show()