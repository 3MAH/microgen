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
supported by the Meshers expression compiler to use compiled evaluation.

.. literalinclude:: ../examples/TPMS/meshers/graded_gyroid.py
   :language: python

Sampled and mapped TPMS
~~~~~~~~~~~~~~~~~~~~~~

The extended backend also accepts offset arrays and ``OffsetGrading`` objects.
Their values are frozen on the original grid, including any normalization.
Tensor grids use trilinear interpolation with linear extrapolation. Clipped
grids use simplex interpolation with nearest-sample extension outside the
sample hull. These paths use NumPy callbacks. ``grading.as_field(grid)``
requires an explicit reference grid for the same reason.

Cylindrical and spherical patches use their original parameter fields and
coordinate maps. Full angular turns are welded at the seam; sector periodicity
exports a rigid ``periodic_transform_*`` matrix instead of a constant shift.
Parameter coordinates and named constraint masks are retained in the VTU file.
Mapped intersections require Tet4 and one thread. On mapped and clipped domains,
Meshers takes ``geometry_tolerance`` in absolute length units by default.
Straight nonperiodic sweeps can adapt their named physical constraints with
Meshers 0.2.1. Spherical domains containing poles or the origin and cylindrical
domains touching their axis are rejected.

``Infill`` and ``GradedInfill`` retain the input object's signed-distance
envelope, which must be closed and manifold. Swept annuli require
``0 < radial_min < radial_max`` and an open regular curve. The forward and
inverse curve maps use the same interpolated frames. Straight Tet4 sweeps use
named-seed adaptation when its options are supported; curved sweeps and other
adaptive options use the Cartesian field route, without mapped periodic ends.
Closed curves and collapsed axes remain unsupported. Quality failures propagate
without changing the requested threshold or falling back to VTK/MMG.

The annular Gyroid regressions cover cell counts ``(16, 8, 48)`` and
``(24, 12, 64)`` with four optimization passes, and scalar ``16`` with one
pass. They check positive tetrahedra, minimum quality 0.1, volume within 1%,
closed surfaces, named constraints and a 250,000-node/tetrahedron limit.
These qualified cases do not establish readiness for arbitrary sweep curves.

.. literalinclude:: ../examples/TPMS/meshers/extended_tpms.py
   :language: python
