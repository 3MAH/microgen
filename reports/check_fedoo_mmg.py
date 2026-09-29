import numpy as np
from types import SimpleNamespace
from convergence import solve

z = np.load("reports/fedoo_debug_mmg_16.npz")
m = SimpleNamespace(
    points=z["points"],
    tetrahedra=z["tetrahedra"],
    surface=z["surface"],
    diagnostics={"minimum_mmg_quality": float(z["quality"].min())},
)
try:
    print("independent", solve(m, affine=True))
except AssertionError as e:
    print("independent affine failure", e)
f = np.sort(
    np.concatenate(
        [m.tetrahedra[:, idx] for idx in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))]
    ),
    axis=1,
)
f, c = np.unique(f, axis=0, return_counts=True)
print(
    "face incidence",
    np.unique(c, return_counts=True),
    "duplicate tets",
    len(m.tetrahedra) - len(np.unique(np.sort(m.tetrahedra, axis=1), axis=0)),
)
p = m.points
t = m.tetrahedra
raw = np.concatenate(
    [t[:, idx] for idx in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))]
)
apices = np.concatenate([t[:, i] for i in (3, 2, 1, 0)])
faces, inverse, count = np.unique(
    np.sort(raw, axis=1), axis=0, return_inverse=True, return_counts=True
)
order = np.argsort(inverse)
ends = np.cumsum(count)
starts = ends - count
sel = np.flatnonzero(count == 2)
a = order[starts[sel]]
b = order[starts[sel] + 1]
f = p[faces[sel]]
normal = np.cross(f[:, 1] - f[:, 0], f[:, 2] - f[:, 0])
d1 = np.einsum("ij,ij->i", normal, p[apices[a]] - f[:, 0])
d2 = np.einsum("ij,ij->i", normal, p[apices[b]] - f[:, 0])
print("same-side adjacent tetrahedra", np.sum(d1 * d2 > 0))
