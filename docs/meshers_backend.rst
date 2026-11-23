.. _meshers_backend:

Volume meshes with Meshers
==========================

Install the optional backend with ``pip install 'microgen[meshers]'``.
These examples export tetrahedral volume meshes as VTU files, retaining
periodic node pairs and boundary facets.

For a ``Tpms`` with a constant offset, the mesh is periodic on the three axes
by default, and ``geometry_tolerance`` is relative to the smallest
``cell_size`` (pass ``geometry_tolerance_mode="absolute"`` for a tolerance in
length units). With ``density``, the offset is fitted on the VTK grid, so the
meshed volume fraction can differ from ``density`` by about 1 %; a larger
``algo_resolution`` reduces the gap.

The uniform Gyroid unit cell is periodic along all three axes:

.. literalinclude:: ../examples/TPMS/meshers/gyroid.py
   :language: python

A pointwise callable can grade the thickness offset. This example varies
the offset along z, so only x and y remain periodic. A linear offset does
not imply a linear density gradient. The callable must use operations
supported by the Meshers expression compiler; sampled arrays and
``OffsetGrading`` objects are not supported by this backend.

.. literalinclude:: ../examples/TPMS/meshers/graded_gyroid.py
   :language: python
