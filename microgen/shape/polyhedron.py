"""Polyhedron.

=============================================
Polyhedron (:mod:`microgen.shape.polyhedron`)
=============================================
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pyvista as pv

from microgen.operations import rotate

from .shape import Shape


if TYPE_CHECKING:
    from microgen.cad import CadShape
    from microgen.shape import KwargsGenerateType, Vector3DType

Vertex = tuple[float, float, float]
Face = dict[str, list[int]]


class Polyhedron(Shape):
    """Class to generate a Polyhedron with a given set of faces and vertices.

    The implicit field is the half-space intersection SDF for a **convex**
    polyhedron: ``f(p) = max_k (n_k · (p - c_k))``, where ``n_k`` is the
    outward unit normal of face ``k`` and ``c_k`` its centroid.  Points
    inside all half-spaces satisfy ``f(p) < 0``.  Voronoi cells produced
    by Neper are convex by construction, which is the primary use case.

    For non-convex polyhedra the implicit field is meaningless; use
    :meth:`generate_cad` / :meth:`generate_surface_mesh` instead (those
    work from the explicit face list).

    .. jupyter-execute::
       :hide-code:

       import microgen

       shape = microgen.Polyhedron().generate_surface_mesh()
       shape.plot(color='white')
    """

    def __init__(
        self: Polyhedron,
        dic: dict[str, list[Vertex | Face]] | None = None,
        **kwargs: Vector3DType,
    ) -> None:
        """Initialize the polyhedron.

        .. warning:
            Give a center parameter only if the polyhedron must be translated
            from its original position.
        """
        super().__init__(**kwargs)
        if dic is None:
            self.dic: dict[str, list[Vertex | Face]] = {
                "vertices": [
                    (1.0, 1.0, 1.0),
                    (1.0, -1.0, -1.0),
                    (-1.0, 1.0, -1.0),
                    (-1.0, -1.0, 1.0),
                ],
                "faces": [
                    {"vertices": [0, 1, 2]},
                    {"vertices": [0, 3, 1]},
                    {"vertices": [0, 2, 3]},
                    {"vertices": [1, 2, 3]},
                ],
            }
        else:
            self.dic = dic

        self.faces_ixs = [face["vertices"] for face in self.dic["faces"]]
        for ixs in self.faces_ixs:
            ixs.append(ixs[0])

        self._setup_frep_field()

    def _setup_frep_field(self: Polyhedron) -> None:
        """Build the convex-polyhedron half-space SDF in world coordinates.

        Bakes the ``center`` translation and ``orientation`` rotation into
        the cached face normals/centroids so the implicit ``_func`` does no
        per-call transform.  Also derives ``_bounds`` from the world-frame
        vertex AABB.
        """
        verts_local = np.asarray(self.dic["vertices"], dtype=np.float64)
        rot = self.orientation.as_matrix()
        cx, cy, cz = (float(c) for c in self.center)
        # World-frame vertices: rotate about origin, then translate to center.
        verts_world = verts_local @ rot.T + np.array([cx, cy, cz])
        # Anchor for outward-orientation flip — the centroid of the convex
        # hull is strictly interior, so any face normal that *should* be
        # outward gives a positive ``n · (face_centroid - poly_centroid)``.
        poly_centroid = verts_world.mean(axis=0)

        normals: list[np.ndarray] = []
        centroids: list[np.ndarray] = []
        for ixs in self.faces_ixs:
            # ``ixs`` is closed (last == first); strip the duplicate.
            face_verts = verts_world[np.asarray(ixs[:-1], dtype=np.int64)]
            if len(face_verts) < 3:
                continue
            edge1 = face_verts[1] - face_verts[0]
            edge2 = face_verts[2] - face_verts[0]
            n = np.cross(edge1, edge2)
            n_len = float(np.linalg.norm(n))
            if n_len == 0.0:
                continue
            n = n / n_len
            face_centroid = face_verts.mean(axis=0)
            # Flip if the cross-product convention pointed inward.
            if np.dot(n, face_centroid - poly_centroid) < 0.0:
                n = -n
            normals.append(n)
            centroids.append(face_centroid)

        if not normals:
            return  # degenerate polyhedron — leave _func unset

        normals_arr = np.asarray(normals, dtype=np.float64)
        centroids_arr = np.asarray(centroids, dtype=np.float64)
        n_dot_c = (normals_arr * centroids_arr).sum(axis=1)

        def _field(
            x: np.ndarray,
            y: np.ndarray,
            z: np.ndarray,
        ) -> np.ndarray:
            x_arr = np.asarray(x)
            shape = x_arr.shape
            p = np.stack(
                [x_arr.ravel(), np.asarray(y).ravel(), np.asarray(z).ravel()],
                axis=-1,
            )
            # signed_d_k(p) = n_k · p − n_k · c_k; shape (N, K)
            sd = p @ normals_arr.T - n_dot_c
            return sd.max(axis=1).reshape(shape)

        self._func = _field
        margin = (
            0.1
            * float(np.linalg.norm(verts_world.max(axis=0) - verts_world.min(axis=0)))
            or 0.1
        )
        vmin = verts_world.min(axis=0) - margin
        vmax = verts_world.max(axis=0) + margin
        self._bounds = (
            float(vmin[0]),
            float(vmax[0]),
            float(vmin[1]),
            float(vmax[1]),
            float(vmin[2]),
            float(vmax[2]),
        )

    def generate_cad(self: Polyhedron, **_: KwargsGenerateType) -> CadShape:
        """Generate a polyhedron CAD shape (OCCT).  Requires the ``[cad]`` extra."""
        from microgen.cad import make_polyhedron

        shape = make_polyhedron(
            vertices=self.dic["vertices"],
            faces_ixs=self.faces_ixs,
            center=self.center,
        )
        return rotate(shape, self.center, self.orientation)

    def generate_surface_mesh(self: Polyhedron, **_: KwargsGenerateType) -> pv.PolyData:
        """Generate a polyhedron VTK shape using the given parameters."""
        faces_pv = copy.deepcopy(self.faces_ixs)
        for vertices_in_face in faces_pv:
            del vertices_in_face[-1]
            vertices_in_face.insert(0, len(vertices_in_face))

        vertices = np.array(self.dic["vertices"])
        faces = np.hstack(faces_pv)

        return pv.PolyData(vertices, faces).compute_normals()


def read_obj(filename: str) -> dict[str, list[Vertex | Face]]:
    """Read vertices and faces from obj format file for polyhedron."""
    dic: dict[str, list[Vertex | Face]] = {"vertices": [], "faces": []}
    with Path(filename).open(encoding="utf-8") as f:
        for line in f:
            data = line.split(" ")
            if data[0] == "v":
                dic["vertices"].append(tuple(float(x) for x in data[1:4]))
            elif data[0] == "f":
                dic["faces"].append(
                    {"vertices": [int(vertex) - 1 for vertex in data[1:]]},
                )
    return dic
