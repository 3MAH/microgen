"""Test utilities for microgen package."""

from typing import Dict, List

import numpy as np


def is_periodic(nodes_coords: np.ndarray, tol: float = 1e-8) -> bool:
    dim = nodes_coords.shape[1]

    min_point = np.min(nodes_coords, axis=0)
    max_point = np.max(nodes_coords, axis=0)

    faces: List[Dict[str, np.ndarray]] = [
        {
            "-": np.where(np.abs(nodes_coords[:, i] - min_point[i]) < tol)[0],
            "+": np.where(np.abs(nodes_coords[:, i] - max_point[i]) < tol)[0],
        }
        for i in range(dim)
    ]

    for i in range(dim):
        if faces[i]["-"].shape != faces[i]["+"].shape:
            return False

    def _sort_dim(
        indices: np.ndarray,
        dim_a: int,
        dim_b: int,
        decimal_round: int = int(-np.log10(tol) - 1),
    ) -> np.ndarray:
        return indices[
            np.lexsort(
                (
                    nodes_coords[indices, dim_a],
                    nodes_coords[indices, dim_b].round(decimal_round),
                )
            )
        ]

    for i in range(dim):
        dim_a = (i + 1) % dim
        dim_b = (i + 2) % dim
        faces[i]["-"] = _sort_dim(faces[i]["-"], dim_a, dim_b)
        faces[i]["+"] = _sort_dim(faces[i]["+"], dim_a, dim_b)

    for i, face in enumerate(faces):
        mask = np.ones(dim, dtype=bool)
        mask[i] = False
        diff = nodes_coords[face["+"]] - nodes_coords[face["-"]]
        if (diff[:, mask] > tol).any():
            return False

    return True
