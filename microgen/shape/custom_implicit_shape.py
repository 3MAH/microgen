"""Custom Implicit Shape.

====================================================================
Custom Implicit Shape (:mod:`microgen.shape.custom_implicit_shape`)
====================================================================

.. jupyter-execute::
   :hide-code:

   pyvista.global_theme.smooth_shading = True
   pyvista.global_theme.split_sharp_edges = True

"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Callable, Literal, Sequence

import cadquery as cq
import numpy as np
import numpy.typing as npt
import pyvista as pv
from scipy.optimize import root_scalar

from .shape import Field, ImplicitShape, Vector3DType


class CustomImplicitShape(ImplicitShape):
    """Custom implicit shape defined by a scalar field.

    :param field: scalar field function that defines the shape
    :param bounds: bounds of the shape
    :param center: center of the shape
    :param orientation: orientation of the shape
    """

    def __init__(
        self: CustomImplicitShape,
        surface_function: Field,
        center: Vector3DType = (0, 0, 0),
        orientation: Vector3DType | None = None,
    ) -> None:
        """Initialize the custom implicit shape."""
        super().__init__(center=center, orientation=orientation)
        self._surface_function = surface_function

    @property
    def surface_function(self: CustomImplicitShape) -> Field:
        """Get the surface function."""
        return self._surface_function
