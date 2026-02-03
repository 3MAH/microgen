.. _RST Tutorials:

Tutorials
=========

This page provides comprehensive tutorials for all Microgen features.

.. jupyter-execute::
   :hide-code:

   import pyvista
   pyvista.set_jupyter_backend('static')
   pyvista.global_theme.background = 'white'
   pyvista.global_theme.window_size = [600, 400]
   pyvista.global_theme.axes.show = False
   pyvista.global_theme.smooth_shading = True
   pyvista.global_theme.split_sharp_edges = True


Basic Shapes
------------

Microgen provides several basic geometric shapes that can be used as building blocks.

Box
^^^

.. jupyter-execute::

   import microgen

   box = microgen.Box(
       center=(0, 0, 0),
       dim=(1.0, 0.6, 0.4),
       orientation=(0, 0, 0)
   )
   box.generate_vtk().plot(color='white')


Sphere
^^^^^^

.. jupyter-execute::

   sphere = microgen.Sphere(
       center=(0, 0, 0),
       radius=0.5
   )
   sphere.generate_vtk().plot(color='white')


Cylinder
^^^^^^^^

.. jupyter-execute::

   cylinder = microgen.Cylinder(
       center=(0, 0, 0),
       height=1.0,
       radius=0.3,
       orientation=(0, 0, 0)
   )
   cylinder.generate_vtk().plot(color='white')


Ellipsoid
^^^^^^^^^

.. jupyter-execute::

   ellipsoid = microgen.Ellipsoid(
       center=(0, 0, 0),
       radii=(0.5, 0.3, 0.2),
       orientation=(0, 0, 0)
   )
   ellipsoid.generate_vtk().plot(color='white')


Capsule
^^^^^^^

.. jupyter-execute::

   capsule = microgen.Capsule(
       center=(0, 0, 0),
       height=1.0,
       radius=0.2,
       orientation=(0, 0, 0)
   )
   capsule.generate_vtk().plot(color='white')


Extruded Polygon
^^^^^^^^^^^^^^^^

.. jupyter-execute::

   # Define polygon corners
   corners = [(0, 0), (1, 0), (1, 1), (0.5, 1.5), (0, 1)]

   polygon = microgen.ExtrudedPolygon(
       center=(0, 0, 0),
       listCorners=corners,
       height=0.5,
       orientation=(0, 0, 0)
   )
   polygon.generate_vtk().plot(color='white')


TPMS Surfaces
-------------

Microgen supports many TPMS (Triply Periodic Minimal Surfaces) functions.

Available Surface Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following TPMS surface functions are available in ``microgen.surface_functions``:

.. list-table:: 3D TPMS Functions
   :widths: 30 70
   :header-rows: 1

   * - Function
     - Description
   * - ``gyroid``
     - Gyroid surface
   * - ``schwarz_p``
     - Schwarz P surface
   * - ``schwarz_d``
     - Schwarz D surface
   * - ``neovius``
     - Neovius surface
   * - ``schoen_iwp``
     - Schoen IWP surface
   * - ``schoen_frd``
     - Schoen FRD surface
   * - ``fischer_koch_s``
     - Fischer-Koch S surface
   * - ``pmy``
     - PMY surface
   * - ``honeycomb``
     - Honeycomb surface
   * - ``lidinoid``
     - Lidinoid surface
   * - ``split_p``
     - Split P surface

.. list-table:: 2D Honeycomb Functions
   :widths: 30 70
   :header-rows: 1

   * - Function
     - Description
   * - ``honeycomb_gyroid``
     - 2D Honeycomb Gyroid
   * - ``honeycomb_schwarz_p``
     - 2D Honeycomb Schwarz P
   * - ``honeycomb_schwarz_d``
     - 2D Honeycomb Schwarz D
   * - ``honeycomb_schoen_iwp``
     - 2D Honeycomb Schoen IWP
   * - ``honeycomb_lidinoid``
     - 2D Honeycomb Lidinoid


TPMS Examples
^^^^^^^^^^^^^

Schwarz P:

.. jupyter-execute::

   schwarz_p = microgen.Tpms(
       surface_function=microgen.surface_functions.schwarz_p,
       offset=0.3,
       cell_size=1.0,
       repeat_cell=2,
       resolution=20
   )
   schwarz_p.sheet.plot(color='white')

Schwarz D:

.. jupyter-execute::

   schwarz_d = microgen.Tpms(
       surface_function=microgen.surface_functions.schwarz_d,
       offset=0.3,
       cell_size=1.0,
       repeat_cell=2,
       resolution=20
   )
   schwarz_d.sheet.plot(color='white')

Neovius:

.. jupyter-execute::

   neovius = microgen.Tpms(
       surface_function=microgen.surface_functions.neovius,
       offset=0.3,
       cell_size=1.0,
       repeat_cell=2,
       resolution=20
   )
   neovius.sheet.plot(color='white')


Using Density Instead of Offset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can specify a target density instead of an offset value:

.. jupyter-execute::

   gyroid_50 = microgen.Tpms(
       surface_function=microgen.surface_functions.gyroid,
       density=0.5,  # 50% density
       cell_size=1.0,
       repeat_cell=2,
       resolution=20
   )
   gyroid_50.sheet.plot(color='white')


Spherical TPMS
^^^^^^^^^^^^^^

Create TPMS on a spherical coordinate system:

.. jupyter-execute::

   spherical_gyroid = microgen.SphericalTpms(
       radius=2.0,
       surface_function=microgen.surface_functions.gyroid,
       offset=0.3,
       cell_size=0.5,
       repeat_cell=(3, 0, 0),  # 0 = auto-fill to complete sphere
       resolution=20
   )
   spherical_gyroid.sheet.plot(color='white')


Cylindrical TPMS
^^^^^^^^^^^^^^^^

Create TPMS on a cylindrical coordinate system:

.. jupyter-execute::

   cylindrical_gyroid = microgen.CylindricalTpms(
       radius=1.5,
       surface_function=microgen.surface_functions.gyroid,
       offset=0.3,
       cell_size=0.5,
       repeat_cell=(2, 0, 3),  # 0 = auto-fill circumference
       resolution=20
   )
   cylindrical_gyroid.sheet.plot(color='white')


3D Operations
-------------

Repeat
^^^^^^

Repeat a shape in a grid pattern:

.. jupyter-execute::

   import microgen
   import pyvista as pv
   import cadquery as cq

   # Create a sphere
   sphere = microgen.Sphere(center=(0, 0, 0), radius=0.15)
   sphere_cad = sphere.generate()

   # Create RVE (Representative Volume Element)
   rve = microgen.Rve(dim=(0.5, 0.5, 0.5))

   # Repeat the sphere in a 3x3x3 grid
   repeated = microgen.repeat_shape(
       unit_geom=sphere_cad,
       rve=rve,
       grid=(3, 3, 3)
   )

   # Export to STL
   cq.exporters.export(repeated, "repeated.stl")

   # Visualize
   plotter = pv.Plotter()
   plotter.add_mesh(pv.read("repeated.stl"))
   plotter.show()


Raster
^^^^^^

Raster a phase to create a periodic pattern:

.. jupyter-execute::

   import cadquery as cq

   # Create ellipsoid and convert to CadQuery shape
   ellipsoid = microgen.Ellipsoid(
       center=(0, 0, 0),
       radii=(0.15, 0.1, 0.08)
   )
   ellipsoid_cad = ellipsoid.generate()
   phase = microgen.Phase(shape=ellipsoid_cad)

   # Define RVE
   rve = microgen.Rve(dim=(0.5, 0.5, 0.5))

   # Raster the phase
   rastered = microgen.raster_phase(
       phase=phase,
       rve=rve,
       grid=(2, 2, 2)
   )

   # Create compound from all solids and export to STL
   compound = cq.Compound.makeCompound(
       [solid for p in rastered for solid in p.solids]
   )
   cq.exporters.export(compound, "rastered.stl")

   # Visualize
   plotter = pv.Plotter()
   plotter.add_mesh(pv.read("rastered.stl"), color='white')
   plotter.show()


Meshing
-------

Regular Mesh with Gmsh
^^^^^^^^^^^^^^^^^^^^^^

Generate a tetrahedral mesh using Gmsh:

.. code-block:: python

   import microgen
   import cadquery as cq

   # Create a TPMS geometry
   gyroid = microgen.Tpms(
       surface_function=microgen.surface_functions.gyroid,
       offset=0.3,
       resolution=20
   )
   shape = gyroid.generate(type_part='sheet')

   # Export to STEP file first
   cq.exporters.export(cq.Workplane().add(shape), 'gyroid.step')

   # Create mesh from STEP file
   microgen.mesh(
       mesh_file='gyroid.step',
       size=0.1,
       order=1,
       output_file='gyroid_mesh.msh'
   )


Periodic Mesh
^^^^^^^^^^^^^

Generate a periodic mesh suitable for homogenization:

.. code-block:: python

   import cadquery as cq

   # Create geometry
   gyroid = microgen.Tpms(
       surface_function=microgen.surface_functions.gyroid,
       offset=0.3,
       cell_size=1.0,
       resolution=20
   )
   shape = gyroid.generate(type_part='sheet')

   # Create a phase from the shape
   phase = microgen.Phase(shape=shape)

   # Export to STEP file first
   cq.exporters.export(cq.Workplane().add(shape), 'gyroid.step')

   # Create periodic mesh
   microgen.mesh_periodic(
       mesh_file='gyroid.step',
       rve=microgen.Rve(dim=(1, 1, 1)),
       list_phases=[phase],
       size=0.05,
       order=1,
       output_file='gyroid_periodic.msh'
   )


Remeshing with Mmg
^^^^^^^^^^^^^^^^^^

Refine or adapt an existing mesh using Mmg:

.. note::

   Mmg must be installed separately to use remeshing features.

.. code-block:: python

   # Use Mmg for remeshing
   mmg = microgen.Mmg(
       input_mesh='input.mesh',
       output_mesh='output.mesh',
       hausd=0.01,  # Hausdorff distance
       hmin=0.01,   # Minimum edge size
       hmax=0.1     # Maximum edge size
   )
   mmg.run()


Further Resources
-----------------

- `Microgen GitHub Repository <https://github.com/3MAH/microgen>`_
- `CadQuery Documentation <https://cadquery.readthedocs.io/en/latest/>`_
- `PyVista Documentation <https://docs.pyvista.org/>`_
- `3MAH Website <https://3mah.github.io/>`_
