.. _RST QuickStart:

Quick-Start Guide
=================

This guide will help you get started with Microgen quickly. We'll cover the basics of creating shapes, TPMS structures, and exporting your results.

.. jupyter-execute::
   :hide-code:

   import pyvista
   pyvista.set_jupyter_backend('static')
   pyvista.global_theme.background = 'white'
   pyvista.global_theme.window_size = [600, 400]
   pyvista.global_theme.axes.show = False
   pyvista.global_theme.smooth_shading = True
   pyvista.global_theme.split_sharp_edges = True


First Steps
-----------

Import microgen:

.. jupyter-execute::

   import microgen


Creating a Basic TPMS Shape
---------------------------

Let's create a Gyroid TPMS structure, one of the most common triply periodic minimal surfaces:

.. jupyter-execute::

   gyroid = microgen.Tpms(
       surface_function=microgen.surface_functions.gyroid,
       offset=0.3,
       cell_size=1.0,
       repeat_cell=2,
       resolution=20
   )

   # Get the sheet geometry as a PyVista mesh
   sheet = gyroid.sheet
   sheet.plot(color='white')


TPMS Part Types
---------------

Each TPMS can generate different part types:

.. jupyter-execute::

   tpms = microgen.Tpms(
       surface_function=microgen.surface_functions.gyroid,
       offset=0.5,
       resolution=20
   )

   # Sheet (wall) geometry
   sheet = tpms.sheet

   # Upper skeletal (network) geometry
   upper_skeletal = tpms.upper_skeletal

   # Lower skeletal (complementary network) geometry
   lower_skeletal = tpms.lower_skeletal

   sheet.plot(color='lightblue')


Creating Basic Shapes
---------------------

Microgen provides several basic geometric shapes:

.. jupyter-execute::

   # Create a sphere
   sphere = microgen.Sphere(center=(0, 0, 0), radius=0.5)
   sphere_vtk = sphere.generate_vtk()
   sphere_vtk.plot(color='white')

.. jupyter-execute::

   # Create a box
   box = microgen.Box(center=(0, 0, 0), dim=(1, 0.5, 0.3))
   box_vtk = box.generate_vtk()
   box_vtk.plot(color='white')

.. jupyter-execute::

   # Create a cylinder
   cylinder = microgen.Cylinder(center=(0, 0, 0), height=1.0, radius=0.3)
   cylinder_vtk = cylinder.generate_vtk()
   cylinder_vtk.plot(color='white')


Generating CAD Geometry
-----------------------

To create CAD-compatible solids (for 3D printing or CAD software):

.. code-block:: python

   # Generate the CadQuery shape
   cad_shape = gyroid.generate(type_part='sheet')


Exporting Results
-----------------

Export to STL (from PyVista mesh):

.. code-block:: python

   sheet.save('gyroid_sheet.stl')

Export to STEP (from CadQuery shape):

.. code-block:: python

   import cadquery as cq

   cad_shape = gyroid.generate(type_part='sheet')
   cq.exporters.export(cq.Workplane().add(cad_shape), 'gyroid_sheet.step')

Export to VTK:

.. code-block:: python

   sheet.save('gyroid_sheet.vtk')


Next Steps
----------

- Explore the :ref:`RST Geometries` for more shape examples
- Check the :ref:`RST Tpms` for all available TPMS surfaces
- Learn about :ref:`RST 3Dop` for repeating and rastering
- See :ref:`RST Mesh` for meshing capabilities
