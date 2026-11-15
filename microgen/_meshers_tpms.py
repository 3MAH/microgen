"""Named domain constraints for mapped and clipped TPMS volume meshes."""

from __future__ import annotations

import numpy as np

from ._meshers import generate_volume
from .shape.tpms import CylindricalTpms, Infill, SphericalTpms, Sweep


def generate_tpms(shape, part, options):
    """Retain the original parameter field and the actual physical envelope."""
    options = dict(options)
    if isinstance(shape, Infill) and not shape.obj.is_manifold:
        raise ValueError(
            "Infill requires a closed manifold envelope for signed distance"
        )
    if isinstance(shape, Sweep):
        if shape.radial_min <= 0:
            raise ValueError(
                "Swept domains require radial_min > 0; the angular map collapses on the axis"
            )
        if np.linalg.norm(shape.curve[0] - shape.curve[-1]) < 1e-9:
            raise NotImplementedError(
                "Closed swept curves require longitudinal seam welding"
            )
    if isinstance(shape, Sweep) and options.get("adapt", False):
        if any(options.get("periodic", (False, False, False))):
            raise NotImplementedError(
                "Adaptive Cartesian sweeps do not support mapped periodic ends"
            )
        raw, envelope = shape.raw_field, shape._cell_box().require_func()
        thickness = shape._offset_func
        scalar = float(shape.offset) if thickness is None else 0.0

        def solid(x, y, z):
            half = (thickness(x, y, z) if thickness is not None else scalar) / 2
            value = raw(x, y, z)
            # Keep both sheet walls distinct during adaptive curvature sizing.
            material = (
                np.maximum(value - half, -value - half)
                if part == "sheet"
                else value + half
                if part == "lower skeletal"
                else half - value
            )
            return np.maximum(material, envelope(x, y, z))

        return generate_volume(
            solid,
            shape.bounds,
            options,
            rotation=shape.orientation,
            center=shape.center,
        )
    mapped = isinstance(shape, (CylindricalTpms, SphericalTpms, Sweep))
    if not mapped and not isinstance(shape, Infill):
        raise NotImplementedError("This TPMS subclass has no Meshers domain adapter")
    raw = shape._parameter_field if mapped else shape.raw_field
    if mapped:
        half = shape.cell_size * shape.repeat_cell / 2
        bounds = tuple(np.column_stack((-half, half)).ravel())
        mapping = shape._map_coordinates
        if shape._parameter_offset_func is not None:
            offset = shape._parameter_offset_func
        elif shape._offset_func is not None:
            physical_offset = shape._offset_func
            physical_map = shape._map_coordinates

            def offset(x, y, z):
                return physical_offset(*physical_map(x, y, z))
        else:

            def offset(x, y, z):
                return float(shape.offset) + 0 * x
    else:
        bounds, mapping = shape.bounds, None
        if shape._offset_func is not None:
            offset = shape._offset_func
        else:

            def offset(x, y, z):
                return float(shape.offset) + 0 * x

    constraints = {}
    if shape.density != 1:
        if part == "sheet":
            constraints["sheet_upper"] = lambda x, y, z: (
                raw(x, y, z) - offset(x, y, z) / 2
            )
            constraints["sheet_lower"] = lambda x, y, z: (
                -raw(x, y, z) - offset(x, y, z) / 2
            )
        elif part == "lower skeletal":
            constraints["lower_skeletal"] = lambda x, y, z: (
                raw(x, y, z) + offset(x, y, z) / 2
            )
        else:
            constraints["upper_skeletal"] = lambda x, y, z: (
                offset(x, y, z) / 2 - raw(x, y, z)
            )
    if not mapped:
        constraints["envelope"] = shape._cell_box().require_func()
    periodic = list(options.get("periodic", (False, False, False)))
    transforms = {}
    weld_axis = None
    if mapped:
        if periodic[0] and not isinstance(shape, Sweep):
            raise ValueError(
                "Radial faces are not related by a rigid periodic transform"
            )
        if isinstance(shape, CylindricalTpms):
            if shape.cylinder_radius <= half[0]:
                raise ValueError(
                    "Cylindrical parameter domains must exclude the singular axis"
                )
            angle = 2 * half[1] / shape.cylinder_radius
            angular_axis = 1
            if periodic[2]:
                transforms[2] = np.eye(4)
                transforms[2][2, 3] = 2 * half[2]
        elif isinstance(shape, SphericalTpms):
            if (
                shape.sphere_radius <= half[0]
                or 2 * half[1] / shape.sphere_radius >= np.pi - 1e-12
            ):
                raise ValueError(
                    "Spherical parameter domains must exclude poles and the origin"
                )
            if periodic[1]:
                raise ValueError(
                    "Polar faces do not have a single rigid periodic transform"
                )
            angle = 2 * half[2] / shape.sphere_radius
            angular_axis = 2
        else:
            if shape.radial_min <= 0:
                raise ValueError(
                    "Swept parameter domains require radial_min > 0; the axis collapses the angular map"
                )
            if periodic[1]:
                raise ValueError("Radial faces are not rigid periodic images")
            if periodic[0]:
                c, n, t = shape._curve_frame_at(
                    np.array([0.0, shape._curve_total_length])
                )
                frames = np.stack((n, np.cross(t, n), t), axis=2)
                transform = np.eye(4)
                transform[:3, :3] = frames[1] @ frames[0].T
                transform[:3, 3] = c[1] - transform[:3, :3] @ c[0]
                transforms[0] = transform
            angle, angular_axis = 2 * np.pi, 2
        if np.isclose(angle, 2 * np.pi, rtol=0, atol=1e-12):
            weld_axis = angular_axis
            periodic[angular_axis] = True
        if periodic[angular_axis]:
            transform = np.eye(4)
            c, s = np.cos(angle), np.sin(angle)
            transform[:2, :2] = [[c, -s], [s, c]]
            transforms[angular_axis] = transform
    options["periodic"] = tuple(periodic)
    if not constraints:
        constraints["domain_lower"] = lambda x, y, z: -1 + 0 * x
        constraints["domain_upper"] = lambda x, y, z: -1 + 0 * x
    elif len(constraints) == 1:
        constraints["domain"] = lambda x, y, z: -1 + 0 * x
    return generate_volume(
        raw,
        bounds,
        options,
        constraints=constraints,
        coordinate_map=mapping,
        periodic_transforms=transforms or None,
        weld_axis=weld_axis,
        rotation=shape.orientation,
        center=shape.center,
    )
