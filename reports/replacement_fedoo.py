"""Fedoo comparison on actual TPMS volume outputs; no shell/solid substitution."""

import json
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
import numpy as np

HERE = Path(__file__).resolve().parent


def child(mode, n):
    import fedoo as fd
    from microgen import Tpms
    from microgen.shape.surface_functions import gyroid
    from microgen.remesh import remesh_keeping_boundaries_for_fem
    from fedoo_convergence import solve

    start = time.perf_counter()
    shape = Tpms(gyroid, offset=0.6, resolution=n)
    if mode == "meshers":
        mesh = shape.generate_meshers(threads=1, periodic=(True,) * 3)
    else:
        grid = remesh_keeping_boundaries_for_fem(
            shape._generate_legacy_volume_mesh(), periodic=True
        )
        points = np.asarray(grid.points, dtype=np.float64).copy()
        tets = grid.cells_dict[10].copy()
        negative = np.linalg.det(points[tets[:, 1:]] - points[tets[:, 0, None]]) < 0
        print(
            json.dumps(
                dict(
                    stage="mmg_mesh",
                    dtype=str(grid.points.dtype),
                    negative_tets=int(negative.sum()),
                    tetrahedra=len(tets),
                )
            ),
            flush=True,
        )
        tets[negative, 0], tets[negative, 1] = (
            tets[negative, 1].copy(),
            tets[negative, 0].copy(),
        )
        faces = np.concatenate(
            [tets[:, idx] for idx in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))]
        )
        faces, counts = np.unique(np.sort(faces, axis=1), axis=0, return_counts=True)
        p = points[tets]
        det = np.linalg.det(p[:, 1:] - p[:, :1])
        edges = sum(
            np.sum((p[:, i] - p[:, j]) ** 2, axis=1)
            for i in range(4)
            for j in range(i + 1, 4)
        )
        quality = np.sqrt(432 * det**2 / edges**3)
        mesh = SimpleNamespace(
            points=points,
            tetrahedra=tets,
            surface=faces[counts == 1],
            diagnostics={"minimum_mmg_quality": float(quality.min())},
        )
        np.savez(
            HERE / f"fedoo_debug_{mode}_{n}.npz",
            points=points,
            tetrahedra=tets,
            surface=mesh.surface,
            quality=quality,
        )
    generation = time.perf_counter() - start
    mms = solve(mesh)
    patch = solve(mesh, affine=True)
    start = time.perf_counter()
    fd.ModelingSpace("3D")
    fem = fd.Mesh(mesh.points.copy(), mesh.tetrahedra.copy(), "tet4")
    assembly = fd.Assembly.create(
        fd.weakform.StressEquilibrium(fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)),
        fem,
        n_elm_gp=1,
    )
    problem = fd.problem.Linear(assembly)
    bottom = np.flatnonzero(np.isclose(mesh.points[:, 2], -0.5, atol=1e-8, rtol=0))
    top = np.flatnonzero(np.isclose(mesh.points[:, 2], 0.5, atol=1e-8, rtol=0))
    assert len(bottom) and len(top)
    problem.bc.add("Dirichlet", bottom, "Disp", 0)
    problem.bc.add("Dirichlet", top, "DispZ", 0.001)
    iterations = [0]

    def count(_):
        iterations[0] += 1

    problem.set_solver(
        "cg", precond=True, rtol=1e-10, atol=1e-12, maxiter=20000, callback=count
    )
    problem.solve()
    u = problem.get_disp().ravel()
    matrix = problem.get_A()
    forces = matrix @ u
    locked = np.concatenate(
        [bottom + i * len(mesh.points) for i in range(3)] + [top + 2 * len(mesh.points)]
    )
    free = np.ones(len(u), bool)
    free[locked] = False
    scale = np.linalg.norm(forces[locked])
    residual = np.linalg.norm(forces[free]) / scale
    assert residual < 1e-7, (iterations, residual)
    energy = float(0.5 * u @ forces)
    reaction = float(forces[top + 2 * len(mesh.points)].sum())
    assert energy > 0 and abs(energy - 0.5 * 0.001 * reaction) / energy < 1e-7
    print(
        json.dumps(
            dict(
                mode=mode,
                resolution=n,
                generation_seconds=generation,
                mms=mms,
                patch=patch,
                compression=dict(
                    stiffness=2 * energy / 0.001**2,
                    energy=energy,
                    reaction=reaction,
                    relative_residual=float(residual),
                    cg_iterations=iterations[0],
                    seconds=time.perf_counter() - start,
                ),
                tetrahedra=len(mesh.tetrahedra),
                points=len(mesh.points),
                minimum_quality=mesh.diagnostics["minimum_mmg_quality"],
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    if len(sys.argv) == 3:
        child(sys.argv[1], int(sys.argv[2]))
    else:
        rows = []
        for mode, n in [
            ("meshers", 12),
            ("mmg", 12),
            ("meshers", 16),
            ("mmg", 16),
            ("meshers", 20),
            ("mmg", 20),
            ("meshers", 28),
        ]:
            p = subprocess.Popen(
                [sys.executable, __file__, mode, str(n)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            try:
                out, err = p.communicate(timeout=180)
                if p.returncode:
                    row = dict(
                        mode=mode,
                        resolution=n,
                        status="error",
                        error=(out + err)[-2500:],
                    )
                else:
                    row = dict(json.loads(out.strip().splitlines()[-1]), status="ok")
            except subprocess.TimeoutExpired:
                subprocess.run(
                    ["taskkill", "/F", "/T", "/PID", str(p.pid)], capture_output=True
                )
                p.communicate()
                row = dict(mode=mode, resolution=n, status="timeout")
            rows.append(row)
            (HERE / "replacement_fedoo.json").write_text(json.dumps(rows, indent=2))
            print(json.dumps(row), flush=True)
