"""Implicit primitive factory functions.

=====================================================================
Implicit Basic Factory (:mod:`microgen.shape.implicit_basic_factory`)
=====================================================================
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .implicit_shape import ImplicitShape


def implicit_sphere(
    center: tuple[float, float, float] = (0, 0, 0),
    radius: float = 1.0,
) -> ImplicitShape:
    """Sphere SDF: ``|p - c| - r``."""
    cx, cy, cz = center
    r = radius
    margin = r * 1.1
    return ImplicitShape(
        func=lambda x, y, z, _cx=cx, _cy=cy, _cz=cz, _r=r: np.sqrt(
            (x - _cx) ** 2 + (y - _cy) ** 2 + (z - _cz) ** 2
        )
        - _r,
        bounds=(
            cx - margin,
            cx + margin,
            cy - margin,
            cy + margin,
            cz - margin,
            cz + margin,
        ),
    )


def implicit_box(
    center: tuple[float, float, float] = (0, 0, 0),
    half_extents: tuple[float, float, float] = (0.5, 0.5, 0.5),
) -> ImplicitShape:
    """Axis-aligned box SDF (exact)."""
    cx, cy, cz = center
    hx, hy, hz = half_extents
    margin = max(hx, hy, hz) * 0.1

    def _box_sdf(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        qx = np.abs(x - cx) - hx
        qy = np.abs(y - cy) - hy
        qz = np.abs(z - cz) - hz
        outside = np.sqrt(
            np.maximum(qx, 0.0) ** 2
            + np.maximum(qy, 0.0) ** 2
            + np.maximum(qz, 0.0) ** 2
        )
        inside = np.minimum(np.maximum(qx, np.maximum(qy, qz)), 0.0)
        return outside + inside

    return ImplicitShape(
        func=_box_sdf,
        bounds=(
            cx - hx - margin,
            cx + hx + margin,
            cy - hy - margin,
            cy + hy + margin,
            cz - hz - margin,
            cz + hz + margin,
        ),
    )


def implicit_cylinder(
    center: tuple[float, float, float] = (0, 0, 0),
    axis: tuple[float, float, float] = (0, 0, 1),
    radius: float = 0.5,
    height: float = 1.0,
) -> ImplicitShape:
    """Capped cylinder SDF along arbitrary axis."""
    c = np.asarray(center, dtype=np.float64)
    a = np.asarray(axis, dtype=np.float64)
    a = a / np.linalg.norm(a)
    half_h = height / 2.0

    def _cyl_sdf(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        px, py, pz = x - c[0], y - c[1], z - c[2]
        proj = px * a[0] + py * a[1] + pz * a[2]
        perp_x = px - proj * a[0]
        perp_y = py - proj * a[1]
        perp_z = pz - proj * a[2]
        d_radial = np.sqrt(perp_x**2 + perp_y**2 + perp_z**2) - radius
        d_axial = np.abs(proj) - half_h
        outside = np.sqrt(
            np.maximum(d_radial, 0.0) ** 2 + np.maximum(d_axial, 0.0) ** 2
        )
        inside = np.minimum(np.maximum(d_radial, d_axial), 0.0)
        return outside + inside

    extent = max(radius, half_h) * 1.1 + np.linalg.norm(c)
    return ImplicitShape(
        func=_cyl_sdf,
        bounds=(
            c[0] - extent,
            c[0] + extent,
            c[1] - extent,
            c[1] + extent,
            c[2] - extent,
            c[2] + extent,
        ),
    )


def implicit_plane(
    point: tuple[float, float, float] = (0, 0, 0),
    normal: tuple[float, float, float] = (0, 0, 1),
) -> ImplicitShape:
    """Half-space SDF: ``dot(p - point, normal)``."""
    pt = np.asarray(point, dtype=np.float64)
    n = np.asarray(normal, dtype=np.float64)
    n = n / np.linalg.norm(n)

    return ImplicitShape(
        func=lambda x, y, z: (x - pt[0]) * n[0]
        + (y - pt[1]) * n[1]
        + (z - pt[2]) * n[2],
        bounds=None,
    )


def implicit_torus(
    center: tuple[float, float, float] = (0, 0, 0),
    major_r: float = 1.0,
    minor_r: float = 0.25,
) -> ImplicitShape:
    """Torus SDF in XZ plane."""
    cx, cy, cz = center

    def _torus_sdf(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        px, py, pz = x - cx, y - cy, z - cz
        q_xz = np.sqrt(px**2 + pz**2) - major_r
        return np.sqrt(q_xz**2 + py**2) - minor_r

    margin = (major_r + minor_r) * 1.1
    return ImplicitShape(
        func=_torus_sdf,
        bounds=(
            cx - margin,
            cx + margin,
            cy - margin,
            cy + margin,
            cz - margin,
            cz + margin,
        ),
    )


def implicit_capsule(
    start: tuple[float, float, float] = (0, 0, -0.5),
    end: tuple[float, float, float] = (0, 0, 0.5),
    radius: float = 0.25,
) -> ImplicitShape:
    """Capsule SDF (line segment distance minus radius)."""
    a = np.asarray(start, dtype=np.float64)
    b = np.asarray(end, dtype=np.float64)
    ab = b - a
    ab_dot = float(np.dot(ab, ab))

    if ab_dot < 1e-30:
        # Degenerate capsule (start == end): fall back to sphere
        return implicit_sphere(center=tuple(a), radius=radius)

    def _capsule_sdf(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        px, py, pz = x - a[0], y - a[1], z - a[2]
        t = np.clip((px * ab[0] + py * ab[1] + pz * ab[2]) / ab_dot, 0.0, 1.0)
        dx = px - t * ab[0]
        dy = py - t * ab[1]
        dz = pz - t * ab[2]
        return np.sqrt(dx**2 + dy**2 + dz**2) - radius

    mid = (a + b) / 2.0
    extent = np.linalg.norm(ab) / 2.0 + radius
    margin = extent * 1.1
    return ImplicitShape(
        func=_capsule_sdf,
        bounds=(
            mid[0] - margin,
            mid[0] + margin,
            mid[1] - margin,
            mid[1] + margin,
            mid[2] - margin,
            mid[2] + margin,
        ),
    )


def implicit_ellipsoid(
    center: tuple[float, float, float] = (0, 0, 0),
    radii: tuple[float, float, float] = (1.0, 0.5, 0.25),
) -> ImplicitShape:
    """Ellipsoid implicit field: ``sqrt((x/rx)^2 + (y/ry)^2 + (z/rz)^2) - 1``.

    This is an approximate SDF but has the correct zero-level set.
    """
    cx, cy, cz = center
    rx, ry, rz = radii
    margin = max(rx, ry, rz) * 1.1

    def _ellipsoid_field(
        x: npt.NDArray[np.float64],
        y: npt.NDArray[np.float64],
        z: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        return (
            np.sqrt(
                ((x - cx) / rx) ** 2
                + ((y - cy) / ry) ** 2
                + ((z - cz) / rz) ** 2
            )
            - 1.0
        )

    return ImplicitShape(
        func=_ellipsoid_field,
        bounds=(
            cx - margin,
            cx + margin,
            cy - margin,
            cy + margin,
            cz - margin,
            cz + margin,
        ),
    )
