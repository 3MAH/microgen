"""Phase class: a collection of OCCT solids belonging to the same phase.

The CAD path goes through OCCT directly via ``OCP`` (installed as the
``[cad]`` extra — ``cadquery-ocp-novtk``); no ``cadquery`` anywhere.  Shapes are
stored as :class:`microgen.cad.CadShape` and solids as raw OCCT
``TopoDS_Solid``.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from .cad import (
    CadShape,
    enumerate_solids,
    make_compound_from_solids,
    make_plane_face,
    split_shape,
    transform_geometry,
    translate_solid,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

    from OCP.TopoDS import TopoDS_Shape, TopoDS_Solid

    from .rve import Rve

    ShapeLike = CadShape | TopoDS_Shape


def _require_ocp() -> None:
    """Raise :class:`ImportError` with an install hint if OCP isn't available."""
    try:
        import OCP  # noqa: F401, PLC0415
    except ImportError as err:
        err_msg = (
            "This Phase operation requires the CAD extra: pip install 'microgen[cad]'"
        )
        raise ImportError(err_msg) from err


def _to_cad_shape(obj: ShapeLike) -> CadShape:
    """Coerce a ``CadShape`` or raw ``TopoDS_Shape`` into a :class:`CadShape`."""
    if isinstance(obj, CadShape):
        return obj
    return CadShape(obj.wrapped if hasattr(obj, "wrapped") else obj)


class Phase:
    """Phase class: a collection of solids with shared material properties.

    Exposes:

    - :attr:`center_of_mass`
    - :attr:`inertia_matrix`
    - :attr:`shape` (a :class:`~microgen.cad.CadShape`)
    - :attr:`solids` (a list of raw OCCT ``TopoDS_Solid``)

    :param shape: a :class:`~microgen.cad.CadShape` or raw ``TopoDS_Shape``
    :param solids: list of raw OCCT solids
    :param center: center
    :param orientation: orientation
    """

    num_instances = 0

    def __init__(
        self: Phase,
        shape: ShapeLike | None = None,
        solids: list[TopoDS_Solid] | None = None,
        center: tuple[float, float, float] | None = None,
        orientation: tuple[float, float, float] | None = None,
    ) -> None:
        """Initialize the phase object."""
        self._shape: CadShape | None = (
            _to_cad_shape(shape) if shape is not None else None
        )
        self._solids: list[TopoDS_Solid] = solids if solids is not None else []
        self.center = center
        self.orientation = orientation

        if shape is None and solids == []:
            warnings.warn("Empty phase", stacklevel=2)

        self.name = f"Phase_{self.num_instances}"

        self._center_of_mass = None
        self._inertia_matrix = None

        Phase.num_instances += 1

    def get_center_of_mass(
        self: Phase,
        *,
        compute: bool = True,
    ) -> npt.NDArray[np.float64]:
        """Return the center of 'mass' of the phase.

        :param compute: if False and center_of_mass already exists, \
            does not compute it (use carefully)
        """
        if isinstance(self._center_of_mass, np.ndarray) and not compute:
            return self._center_of_mass
        self._center_of_mass = self._compute_center_of_mass()
        return self._center_of_mass

    center_of_mass = property(get_center_of_mass)

    def _compute_center_of_mass(self: Phase) -> npt.NDArray[np.float64]:
        """Calculate the center of 'mass' of the phase."""
        _require_ocp()
        from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
        from OCP.GProp import GProp_GProps  # noqa: PLC0415

        if self.shape is None:
            err_msg = "Cannot compute center of mass on an empty phase"
            raise ValueError(err_msg)

        properties = GProp_GProps()
        BRepGProp.VolumeProperties_s(self.shape.wrapped, properties)

        com = properties.CentreOfMass()
        return np.array([com.X(), com.Y(), com.Z()])

    def get_inertia_matrix(
        self: Phase,
        *,
        compute: bool = True,
    ) -> npt.NDArray[np.float64]:
        """Calculate the inertia matrix of the phase.

        :param compute: if False and inertia_matrix already exists, \
            does not compute it (use carefully)
        """
        if isinstance(self._inertia_matrix, np.ndarray) and not compute:
            return self._inertia_matrix
        self._inertia_matrix = self._compute_inertia_matrix()
        return self._inertia_matrix

    inertia_matrix = property(get_inertia_matrix)

    def _compute_inertia_matrix(self: Phase) -> npt.NDArray[np.float64]:
        """Calculate the inertia matrix of the phase."""
        _require_ocp()
        from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
        from OCP.GProp import GProp_GProps  # noqa: PLC0415

        if self.shape is None:
            err_msg = "Cannot compute inertia matrix on an empty phase"
            raise ValueError(err_msg)

        properties = GProp_GProps()
        BRepGProp.VolumeProperties_s(self.shape.wrapped, properties)

        inm = properties.MatrixOfInertia()
        return np.array(
            [
                [inm.Value(1, 1), inm.Value(1, 2), inm.Value(1, 3)],
                [inm.Value(2, 1), inm.Value(2, 2), inm.Value(2, 3)],
                [inm.Value(3, 1), inm.Value(3, 2), inm.Value(3, 3)],
            ],
        )

    @property
    def shape(self: Phase) -> CadShape | None:
        """Return the shape of the phase as a :class:`~microgen.cad.CadShape`."""
        if self._shape is not None:
            return self._shape
        if len(self._solids) > 0:
            self._shape = make_compound_from_solids(self._solids)
            return self._shape

        warnings.warn("No shape or solids", stacklevel=2)
        return None

    @property
    def solids(self: Phase) -> list[TopoDS_Solid]:
        """Return the list of OCCT ``TopoDS_Solid`` in the phase."""
        if len(self._solids) > 0:
            return self._solids
        if self._shape is not None:
            self._solids = enumerate_solids(self._shape)
            return self._solids

        warnings.warn("No solids or shape", stacklevel=2)
        return []

    def translate(self: Phase, vec: Sequence[float]) -> None:
        """Translate phase by a given vector (in place)."""
        if self._shape is None:
            err_msg = "Cannot translate a phase with no shape"
            raise ValueError(err_msg)
        self._shape = self._shape.translate(vec)
        self._center_of_mass = self._compute_center_of_mass()

    @staticmethod
    def rescale_shape(
        shape: ShapeLike,
        scale: float | tuple[float, float, float],
    ) -> CadShape:
        """Rescale ``shape`` by ``scale`` = ``(sx, sy, sz)`` (or a scalar).

        Preserves the shape's center of mass — scaling is performed about it.
        """
        shape = _to_cad_shape(shape)
        if isinstance(scale, float):
            scale = (scale, scale, scale)

        center = shape.center()
        cx, cy, cz = center.x, center.y, center.z
        sx, sy, sz = (float(s) for s in scale)

        # Equivalent to: translate(-c) → scale about origin → translate(+c)
        # Expressed as a single 3x4 affine matrix for BRepBuilderAPI_GTransform.
        matrix = np.array(
            [
                [sx, 0.0, 0.0, cx - sx * cx],
                [0.0, sy, 0.0, cy - sy * cy],
                [0.0, 0.0, sz, cz - sz * cz],
            ],
            dtype=np.float64,
        )
        return transform_geometry(shape, matrix)

    def rescale(self: Phase, scale: float | tuple[float, float, float]) -> None:
        """Rescale phase (in place) by ``scale = (sx, sy, sz)`` or a scalar."""
        if self._shape is None:
            err_msg = "Cannot rescale a phase with no shape"
            raise ValueError(err_msg)
        self._shape = self.rescale_shape(self._shape, scale)

    @staticmethod
    def repeat_shape(
        unit_geom: ShapeLike,
        rve: Rve,
        grid: tuple[int, int, int],
    ) -> CadShape:
        """Repeat ``unit_geom`` on a ``grid`` within the ``rve`` periodicity cell.

        Returns a :class:`~microgen.cad.CadShape` wrapping an OCCT compound
        of translated copies of ``unit_geom``.
        """
        unit_geom = _to_cad_shape(unit_geom)
        center = np.array(unit_geom.center().to_tuple())

        copies: list[TopoDS_Solid] = []
        for idx in np.ndindex(*grid):
            pos = center - rve.dim * (0.5 * np.array(grid) - 0.5 - np.array(idx))
            copies.append(translate_solid(unit_geom.wrapped, pos))
        return make_compound_from_solids(copies)

    def repeat(self: Phase, rve: Rve, grid: tuple[int, int, int]) -> None:
        """Repeat phase in place on a ``grid`` within the ``rve`` periodicity cell."""
        if self.shape is None:
            err_msg = "Cannot repeat a phase with no shape"
            raise ValueError(err_msg)
        self._shape = self.repeat_shape(self.shape, rve, grid)

    def split_solids(self: Phase, rve: Rve, grid: list[int]) -> list[TopoDS_Solid]:
        """Split solids from phase according to the rve divided by the given grid.

        Each solid is split by the (grid-1) interior planes along each axis;
        planes are constructed with :func:`microgen.cad.make_plane_face` and
        applied through OCCT's ``BRepAlgoAPI_Splitter``.

        :param rve: RVE divided by the given grid
        :param grid: number of divisions in each direction ``[x, y, z]``

        :return: list of raw OCCT ``TopoDS_Solid``
        """
        result: list[TopoDS_Solid] = []
        for solid in self.solids:
            current = CadShape(solid)
            for dim in range(3):
                direction = tuple(int(dim == i) for i in range(3))
                coords = np.linspace(
                    start=rve.min_point[dim],
                    stop=rve.max_point[dim],
                    num=grid[dim],
                    endpoint=False,
                )[1:]
                for pos in coords:
                    base_pnt = tuple(float(pos) * direction[k] for k in range(3))
                    plane = make_plane_face(base_pnt, direction)
                    current = split_shape(current, plane)
            result.extend(enumerate_solids(current))
        return result

    def rasterize(
        self: Phase,
        rve: Rve,
        grid: list[int],
        *,
        phase_per_raster: bool = True,
    ) -> list[Phase] | None:
        """Raster solids from phase according to the rve divided by the given grid.

        :param rve: RVE divided by the given grid
        :param grid: number of divisions in each direction [x, y, z]
        :param phase_per_raster: if True, returns list of phases

        :return: list of Phases if required
        """
        solids = self.split_solids(rve, grid)

        if phase_per_raster:
            return self.generate_phase_per_raster(solids, rve, grid)
        self._solids = solids
        self._shape = make_compound_from_solids(self._solids)
        return None

    @classmethod
    def generate_phase_per_raster(
        cls: type[Phase],
        solids: list[TopoDS_Solid],
        rve: Rve,
        grid: list[int],
    ) -> list[Phase]:
        """Raster solids into per-grid-cell Phases.

        :param solids: list of OCCT solids
        :param rve: RVE divided by the given grid
        :param grid: number of divisions in each direction ``[x, y, z]``

        :return: list of :class:`Phase`
        """
        _require_ocp()
        from OCP.BRepGProp import BRepGProp  # noqa: PLC0415
        from OCP.GProp import GProp_GProps  # noqa: PLC0415

        grid_arr = np.array(grid)
        solids_phases: list[list[TopoDS_Solid]] = [
            [] for _ in range(int(np.prod(grid_arr)))
        ]
        for solid in solids:
            props = GProp_GProps()
            BRepGProp.VolumeProperties_s(solid, props)
            com = props.CentreOfMass()
            center = np.array([com.X(), com.Y(), com.Z()])
            i, j, k = np.floor(
                grid_arr * (center - rve.min_point) / rve.dim,
            ).astype(int)
            ind = i + grid_arr[0] * j + grid_arr[0] * grid_arr[1] * k
            solids_phases[int(ind)].append(solid)
        return [Phase(solids=solids) for solids in solids_phases if len(solids) > 0]
