# Optimizer and GIL experiments

These experiments follow `d4f276f`. They are opt-in; the browser package and normal meshing defaults are unchanged. Raw measurements and output fingerprints are in `optimizer-experiments-results.json`.

## Optimizer variants

Analytical mode differentiates the tetrahedron penalty exactly, including simultaneous translation of periodic copies. It projects that gradient onto the allowed tangent directions. The implicit surface-error term retains its projected finite differences. All candidate acceptance and final geometry checks remain enabled. This changes the proposed direction slightly, so identical final coordinates are not expected.

Active-set mode skips a group only after a zero accepted move and while no incident tetrahedron's vertices have changed. It invalidates neighboring periodic groups after accepted moves and resets at the pass where the step size changes. It currently applies only to legacy serial ordering, `threads=0`, and adds adjacency storage. It is not implemented for the colored parallel path.

Native release builds, three runs per case, median seconds:

| Case | Baseline | Analytical | Active set | Both |
|---|---:|---:|---:|---:|
| Gyroid, one cell | 0.616 | 0.540 | 0.533 | 0.582 |
| Split-P, one cell | 1.642 | 1.574 | 1.493 | 1.583 |
| Schwarz-P, one cell | 0.507 | 0.452 | 0.455 | 0.430 |
| Graded gyroid, 2³ cells | 3.583 | 3.543 | 3.437 | 3.307 |

Unit cells use 16 samples per cell; graded cases use 12. All use offset 0.6, four passes before and after refinement, and tolerance 0.01. Grading is 30% along X. Timings cover generation and adapter metrics but exclude final JSON serialization. Variants ran sequentially on a shared desktop, not interleaved. The initial baseline ran slower; the table uses the complete repeated baseline. Differences of a few percent are inconclusive, and combined changes did not consistently outperform either alone.

The active-set outputs match the baseline fingerprints exactly. Combined outputs match analytical-only fingerprints exactly. All four successful cases retain element counts, closed edges, periodic face matching and FEA targets. Analytical minimum qualities change as follows:

| Case | Baseline minimum quality | Analytical minimum quality |
|---|---:|---:|
| Gyroid | 0.143534 | 0.143536 |
| Split-P | 0.111999 | 0.111999 |
| Schwarz-P | 0.451992 | 0.451994 |
| Graded gyroid | 0.141570 | 0.141519 |

The coarse graded Split-P 2³ case fails in every variant: sampled surface error is about 0.02163 against a 0.01 tolerance. No partial mesh is accepted. The initial baseline runner stopped on this error; subsequent runs record it explicitly. These optimizations do not fix this existing limitation.

Decision: retain both as experiments. Active-set correctness looks promising, but the speedup is modest. Analytical derivatives pass an independent finite-difference test, but keeping the surface penalty numerical limits the work removed. Neither result justifies a broad default change from this sample alone.

## GIL release and parallelism

The Python extension already calls `py.detach` around native generation. Compiled evaluators run without Python callbacks. Callback evaluators reacquire Python through `Python::attach`. Releasing the GIL enables concurrency; it does not make one serial computation faster by itself.

Four independent compiled-field gyroid meshes, each with one native worker, 15 lattice divisions, four optimization passes before/after refinement:

| Python job threads | Total time for four meshes | Relative throughput |
|---|---:|---:|
| 1 | 3.725 s | 1.00× |
| 2 | 1.868 s | 1.99× |
| 4 | 0.920 s | 4.05× |

All jobs produce identical point/tetrahedron fingerprints, 15,909 tetrahedra and minimum quality 0.144037. The 4.05× ratio includes measurement noise; it is approximately fourfold throughput, not a fourfold reduction in individual mesh latency. Compilation is prepared once before timing.

The smaller Python-callback probe uses 10 divisions and one pass before/after refinement. Four jobs take 2.219 / 2.204 / 2.286 seconds at 1 / 2 / 4 Python threads. Each job calls Python 49,252 times, and all output fingerprints match. This probe demonstrates poor callback concurrency; its timings must not be compared directly with the higher-resolution compiled probe.

For one compiled graded 2³ gyroid, 23 divisions and four passes before/after refinement:

| Native workers | Single-mesh time |
|---|---:|
| 1 | 4.472 s |
| 2 | 3.802 s |
| 4 | 5.823 s |

These worker counts all produce identical fingerprints, 63,682 tetrahedra and minimum quality 0.148232. They use the existing colored update order, which differs from the `threads=0` optimizer experiment. Four workers regress on this case. The measurements establish the regression, not its exact breakdown into scheduling, allocation or cache costs.

For microgen, use compiled fields and gradients, reuse prepared evaluators, and parallelize independent jobs with a bounded thread pool and one native worker per job. For one large graded domain, benchmark a small native worker count instead of selecting all cores automatically. Do not multiply a large outer job pool by a large inner worker count. Python callbacks remain useful for unsupported field expressions, but do not assume they scale like compiled fields.

## Reproduction and validation

The Cargo feature `optimizer-experiments` gates the variants. Set `MESHER_OPT_EXPERIMENT` before the process first enters optimization: `baseline`, `analytic`, `active`, or `analytic-active`. The selection is cached for that process. The Python crate exposes the same opt-in Cargo feature; regular wheels do not enable it.

```powershell
$env:MESHER_OPT_EXPERIMENT='analytic-active'
cargo test -p meshers-wasm --release --features meshers-core/optimizer-experiments profile_optimizer_matrix -- --ignored --nocapture
```

Set `MESHER_CPU_PROFILE=1` for optimizer phase timings and active-set visit/skip counts. The GIL benchmark runs against the rebuilt default Python extension, with experiments disabled:

```powershell
python crates/meshers-python/examples/gil_benchmark.py --output gil-results.json
```

The examples assume the experimental source package is on `PYTHONPATH`. Each GIL configuration runs three times and checks output fingerprints. The optimizer comparison script checks exact active-set equivalence and analytical quality/periodicity. The core suite passes 23 tests in both baseline and combined experiment modes, including derivative and active-set invalidation tests.
