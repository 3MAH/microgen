"""Bounded parallel generation of independent native TPMS meshes."""

from concurrent.futures import ThreadPoolExecutor
from operator import index

import meshers


def _available_threads():
    # Use the same affinity-aware, capped CPU budget as meshers threads=None.
    return meshers._meshers.available_threads()


def _integer(value, name, minimum):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    try:
        value = index(value)
    except TypeError as error:
        raise ValueError(f"{name} must be an integer >= {minimum}") from error
    if value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def generate_meshers_parallel(shapes, *, max_workers=None, threads=1, **options):
    """Mesh independent TPMS objects concurrently, preserving input order.

    Uses all available job slots by default, with one native worker per job.
    Compiled fields release the GIL; Python callbacks may limit scaling.
    ``max_workers`` limits simultaneous jobs. ``threads`` limits native workers
    per job; None divides the CPU budget between jobs, and 0 selects meshers'
    serial ordering. The product of job slots and native workers never exceeds
    meshers' available CPU budget for this batch.

    Pass the same meshing options accepted by ``Tpms.generate_meshers``.
    Returns native meshes, not PyVista objects. Each shape must be a distinct
    instance: density fitting can change its offset. Do not mutate shapes or
    shared field state during generation. A failure raises; no partial result
    list is returned. Already-running jobs finish before the call returns.
    Separate simultaneous batch calls have separate budgets.
    """
    items = list(shapes)
    if len({id(shape) for shape in items}) != len(items):
        raise ValueError("Use distinct shape instances for parallel meshing")
    budget = max(1, _available_threads())
    limit = budget if max_workers is None else _integer(max_workers, "max_workers", 1)
    if threads is not None:
        threads = _integer(threads, "threads", 0)
        if threads > budget:
            raise ValueError(f"threads exceeds the available CPU budget ({budget})")
    if not items:
        return []
    workers = min(len(items), limit, budget // max(1, threads or 1))
    native_threads = max(1, budget // workers) if threads is None else threads

    def generate(shape):
        result = shape.generate_meshers(threads=native_threads, **options)
        result.diagnostics["batch_workers"] = workers
        result.diagnostics["batch_native_threads"] = native_threads
        result.diagnostics["batch_cpu_budget"] = budget
        return result

    with ThreadPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(generate, items))
