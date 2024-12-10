"""Grading functions for TPMS."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import pyvista as pv
from typing_extensions import override

if TYPE_CHECKING:
    import numpy as np
    import numpy.typing as npt
    import pyvista as pv


class OffsetGrading(ABC):
    """Base class for offset grading functions."""

    @abstractmethod
    def compute_offset(
        self: OffsetGrading,
        grid: pv.UnstructuredGrid | pv.StructuredGrid,
    ) -> npt.NDArray[np.float64]:
        """Compute the offset of the grid.

        This method should compute the offset on each point of the grid and return \
            it as a 1D array.
        The `lower_surface` and `upper_surface` fields of the grid will then be \
            updated using this computed offset.

        Parameters
        ----------
        grid : pv.UnstructuredGrid | pv.StructuredGrid
            The grid to compute the offset on.

        Returns
        -------
        npt.NDArray[np.float64]
            The offset of the grid as a 1D array that matches the \
                number of points in the grid.

        """
        raise NotImplementedError


class NormedDistance(OffsetGrading):
    """Compute the offset based on the implicit distance to an object."""

    def __init__(
        self: NormedDistance,
        obj: pv.PolyData,
        boundary_offset: float,
        furthest_offset: float,
        boundary_weight: float = 1.0,
    ) -> None:
        """Initialize the ImplicitDistance object.

        Parameters
        ----------
        obj : pv.PolyData
            The object to compute the implicit distance to.
        boundary_offset : float
            The offset applied on the boundary of the object.
        furthest_offset : float
            The offset applied at the furthest point from the object.
        boundary_weight : float, optional
            Coefficient to influence the offset on the boundary. The default is 1.0.
            The normed distance is 0 at the object surface and 1 at the
            furthest point from the object surface.
            A boundary weight of 2.0 would make the boundary thicker than the default.
            A boundary weight of 0.5 would make the boundary thinner than the default.

        """
        self.obj = obj
        self.boundary_offset = boundary_offset
        self.furthest_offset = furthest_offset
        self.boundary_weight = boundary_weight

    @override
    def compute_offset(
        self: NormedDistance,
        grid: pv.UnstructuredGrid | pv.StructuredGrid,
    ) -> npt.NDArray[np.float64]:
        distance = grid.compute_implicit_distance(self.obj)["implicit_distance"]
        max_distance = distance.max()

        if max_distance == 0:
            err_msg = "The maximum distance to the object is 0. \
                Cannot compute the normed distance offset."
            raise ValueError(err_msg)

        normed_distance = distance / max_distance
        normed_distance[normed_distance < 0] = 0

        weighted_distance = normed_distance**self.boundary_weight

        return (
            self.furthest_offset
            + (self.boundary_offset - self.furthest_offset) * weighted_distance
        )
