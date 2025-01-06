"""Tests for the external module."""

import pytest

import microgen
from microgen.external import MmgError


def test_run_mmg2d_with_no_arguments_must_raise_mmg_error() -> None:
    """Test that mmg2d command without arguments raises MmgError."""
    with pytest.raises(MmgError):
        microgen.external.Mmg.mmg2d()


def test_run_mmgs_with_no_arguments_must_raise_mmg_error() -> None:
    """Test that mmgs command without arguments raises MmgError."""
    with pytest.raises(MmgError):
        microgen.external.Mmg.mmgs()


def test_run_mmg3d_with_no_arguments_must_raise_mmg_error() -> None:
    """Test that mmg3d command without arguments raises MmgError."""
    with pytest.raises(MmgError):
        microgen.external.Mmg.mmg3d()
