"""Test the implicit shapes."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

import microgen
from microgen.shape import Box, ImplicitShape


def test_explicit_box_is_implicit_must_return_false() -> None:
    """Test if the explicit box is not an implicit shape."""
    box = Box(dim=(1, 1, 1), implicit=False)
    assert not isinstance(box, ImplicitShape)


def test_implicit_box_is_implicit_must_return_true() -> None:
    """Test if the implicit box is an implicit shape."""
    box = Box(dim=(1, 1, 1), implicit=True)
    assert isinstance(box, ImplicitShape)
