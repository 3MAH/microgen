from itertools import product
from pathlib import Path

import cadquery as cq
import numpy as np

from microgen import CustomLattice

# Generate custom auxetic strut lattice based on
# DOI:10.1007/s00707-019-02387-x and
# https://doi.org/10.1016/j.ijmecsci.2017.05.048

theta = np.pi / 3
l2 = 1.0 / (2.0 * np.sin(theta))  # l
l1 = 0.5 + l2 * np.cos(theta)  # h
l3 = (l1 / 2.0) - np.cos(theta)

outer_cube_vertices = np.array(list(product([-0.5, 0.5], repeat=3)))
z_axis_faces_centers = np.array([[0, 0, z] for z in [0.5, -0.5]])
z_axis_reentrant_vertices = np.array([[0, 0, z] for z in [l3, -l3]])
xy_face_reentrant_vertices = np.array(
    [[x, y, z] for x, y in product([0.5, -0.5], [0.0]) for z in [l1 / 2, -l1 / 2]]
    + [[x, y, z] for x, y in product([0.0], [0.5, -0.5]) for z in [l1 / 2, -l1 / 2]]
)
outer_cube_edges_reentrant_vertices = np.array(
    list(product([0.5, -0.5], [0.5, -0.5], [l3, -l3]))
)
base_vertices = np.vstack(
    [
        outer_cube_vertices,
        outer_cube_edges_reentrant_vertices,
        xy_face_reentrant_vertices,
        z_axis_faces_centers,
        z_axis_reentrant_vertices,
    ]
)

strut_heights = np.array([l1] * 4 + [l2] * 24 + [l3 + 0.5] * 10)

l1_strut_vertex_pairs = np.array([[16, 17], [20, 21], [18, 19], [22, 23]])
l2_strut_vertex_pairs = np.array(
    [
        [18, 15],
        [19, 14],
        [18, 13],
        [12, 19],
        [15, 22],
        [14, 23],
        [23, 10],
        [11, 22],
        [11, 16],
        [10, 17],
        [16, 9],
        [8, 17],
        [9, 20],
        [8, 21],
        [13, 20],
        [12, 21],
        [20, 27],
        [18, 27],
        [27, 22],
        [16, 27],
        [23, 26],
        [17, 26],
        [26, 21],
        [26, 19],
    ]
)
half_l1_strut_vertex_pairs = np.array(
    [
        [9, 7],
        [11, 5],
        [8, 6],
        [10, 4],
        [26, 25],
        [12, 2],
        [13, 3],
        [27, 24],
        [15, 1],
        [14, 0],
    ]
)

strut_vertex_pairs = np.vstack(
    [l1_strut_vertex_pairs, l2_strut_vertex_pairs, half_l1_strut_vertex_pairs]
)

auxetic_lattice = CustomLattice(
    strut_heights=strut_heights,
    base_vertices=base_vertices,
    strut_vertex_pairs=strut_vertex_pairs,
    strut_joints=True,
)

shape = auxetic_lattice.generate()
stl_file = Path(__file__).parent / "auxetic_custom_lattice.stl"
cq.exporters.export(shape, stl_file.name)
