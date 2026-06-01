"""Tests for the ``Phase.pieces`` invariant — the "collection of sub-pieces".

A :class:`Phase` is conceptually "one labelled material region", but
practically it may split into multiple sub-pieces under periodicity or
rasterisation.  :attr:`Phase.pieces` surfaces that structure uniformly
across backends:

- field-backed: ``scipy.ndimage.label`` on ``{field < iso}``.
- CAD-backed: one :class:`Piece` per ``TopoDS_Solid``.
"""

import numpy as np

from microgen import Phase, Sphere
from microgen.phase import Piece
from microgen.shape.implicit_ops import union

# ruff: noqa: S101


def test_pieces_field_single_sphere() -> None:
    """A single sphere has exactly one piece."""
    phase = Phase.from_shape(Sphere(radius=0.4))
    pieces = phase.pieces
    assert len(pieces) == 1
    p = pieces[0]
    assert isinstance(p, Piece)
    assert p.voxel_mask is not None
    assert p.cad is None
    # The volume estimate is correct up to grid resolution.
    expected = (4.0 / 3.0) * np.pi * 0.4**3
    assert np.isclose(p.volume, expected, rtol=0.05)


def test_pieces_field_two_disjoint_spheres() -> None:
    """Two disjoint spheres produce two connected components.

    Built by SDF union — the union field is negative inside *either*
    sphere; the two solid regions are spatially disjoint, so
    ``scipy.ndimage.label`` yields two components.
    """
    a = Sphere(center=(-0.6, 0.0, 0.0), radius=0.2)
    b = Sphere(center=(0.6, 0.0, 0.0), radius=0.2)
    merged = union(a, b)
    # The unioned shape's auto-bounds is the AABB of the two — explicitly
    # bound the phase to span both spheres.
    phase = Phase.from_shape(merged, bounds=(-1.0, 1.0, -0.4, 0.4, -0.4, 0.4))
    pieces = phase.pieces
    assert len(pieces) == 2
    coms = sorted(p.com[0] for p in pieces)
    assert np.isclose(coms[0], -0.6, atol=0.05)
    assert np.isclose(coms[1], 0.6, atol=0.05)


def test_pieces_cad_single_sphere() -> None:
    """A CAD-backed sphere has one piece carrying a CadShape payload."""
    cad = Sphere(radius=0.5).generate_cad()
    phase = Phase.from_cad(cad)
    pieces = phase.pieces
    assert len(pieces) == 1
    assert pieces[0].cad is not None
    assert pieces[0].voxel_mask is None


def test_pieces_cached() -> None:
    """``phase.pieces`` is a ``@cached_property`` — second access is the same list."""
    phase = Phase.from_shape(Sphere(radius=0.4))
    assert phase.pieces is phase.pieces


def test_phase_translated_preserves_period_and_invalidates_cache() -> None:
    """``translated()`` returns a new Phase — pieces and COM are recomputed."""
    a = Sphere(center=(-0.6, 0.0, 0.0), radius=0.2)
    b = Sphere(center=(0.6, 0.0, 0.0), radius=0.2)
    merged = union(a, b)
    phase = Phase.from_shape(merged, bounds=(-1.0, 1.0, -0.4, 0.4, -0.4, 0.4))
    moved = phase.translated((1.0, 0.0, 0.0))
    # New Phase => caches independent; both still report 2 pieces.
    assert len(moved.pieces) == 2
    shifts = sorted(p.com[0] for p in moved.pieces)
    assert np.isclose(shifts[0], -0.6 + 1.0, atol=0.05)
    assert np.isclose(shifts[1], 0.6 + 1.0, atol=0.05)
