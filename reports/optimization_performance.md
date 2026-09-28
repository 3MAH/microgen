# Optimization experiment

Baseline: `345e872`. Measurements are local release builds on the same Windows machine. Each case uses one warm-up and three measured runs; tables show medians. Variants ran sequentially, so small timing differences can reflect machine load. Raw results, including the intermediate analytical-gradient variant, are in `performance-results.json`.

## FEA volume generation

The grid, eight optimization passes, periodic constraints, geometry tolerance and quality targets are unchanged. Analytical TPMS derivatives replace six finite-difference field evaluations per gradient in the browser adapter. The core optimizer stops evaluating a candidate once its tetrahedron quality violates the existing acceptance floor, or its accumulated nonnegative penalty cannot beat the current best cost. Finite-difference cost-gradient probes still receive full evaluations.

| Native, serial | Before | After | Speedup |
|---|---:|---:|---:|
| Gyroid, one cell, 15,909 tetrahedra | 0.902 s | 0.649 s | 1.39× |
| Split-P, one cell, 21,526 tetrahedra | 3.360 s | 1.660 s | 2.02× |

| WASM through the JS adapter | Before | After | Speedup |
|---|---:|---:|---:|
| Gyroid, one cell, 16 samples | 2.120 s | 1.371 s | 1.55× |
| Split-P, one cell, 16 samples | 5.213 s | 2.616 s | 1.99× |
| Graded gyroid, 2³ cells, 12 samples per cell | 10.684 s | 6.606 s | 1.62× |

Adapter timings include meshing, diagnostics and JSON conversion, but exclude module loading, rendering and worker transfer. Native timings include Rust diagnostics and serialization. WASM and native are separate runtime measurements, not interchangeable predictions for every browser.

The three FEA cases retain their element counts, closure and periodic face matches. Minimum MMG quality remains 0.143534, 0.111999 and 0.141570 respectively. Minimum quality, sampled surface error and volume differ by less than 1e-8 from baseline. This is a speed improvement without reducing optimization passes. The remaining cost is still substantial; no sub-100-ms FEA result is claimed.

## Fast triangle extraction

The direct linear extractor now packs edge-cache keys into a single integer and skips cubes whose sampled range cannot produce a wall or cap. Cache memory still scales with crossing edges; this does not allocate a dense edge table. A local diagonal choice maximizes the smallest triangle angle using squared lengths and cross products, without extra field evaluations or a smoothing pass.

The native extraction-only benchmark uses graded geometry along X, 16 samples per cell and no periodic tiling. Larger cases evaluate the full domain. These times exclude diagnostics and serialization.

| Geometry | Cells | Before | After | Triangles |
|---|---:|---:|---:|---:|
| Gyroid | 1³ | 2.291 ms | 2.142 ms | 13,504 |
| Gyroid | 3³ | 66.187 ms | 49.691 ms | 372,744 |
| Gyroid | 7³ | 1.164 s | 1.007 s | 4,755,368 |
| Split-P | 1³ | 4.973 ms | 4.181 ms | 22,228 |
| Split-P | 3³ | 131.443 ms | 103.078 ms | 607,728 |
| Split-P | 7³ | 2.426 s | 1.805 s | 7,756,296 |

The 7³ cases are 14% and 26% faster. The final native surface runs include coarse profiling. VTK was not rerun in this experiment; these are improvements over the previously benchmarked direct extractor, not new VTK speed ratios.

The separate WASM validation cases use uniform one-cell meshes and graded 3³ meshes at 20 samples per cell. All retain counts, closure and periodic matches. The worst and first-percentile triangle angles are unchanged. More triangles enter the 55–60° minimum-angle bin, while area coefficient of variation increases slightly: gyroid 0.8273→0.8288 and Split-P 0.8180→0.8193 in the unit-cell cases. The diagonal choice therefore offers a modest angle-distribution improvement, not sliver removal or FEA-quality surfaces.

## Reproduction

```powershell
cargo test -p meshers-wasm --release profile_ -- --ignored --nocapture --test-threads=1
# Set MESHER_CPU_PROFILE=1 for coarse native phase timings.
cd demos/tpms-browser
node performance.mjs path/to/meshers.wasm results.json
```

Validation covers analytical derivatives against independent central differences, rejection of candidates without unnecessary surface evaluations, core geometry/topology tests, all three TPMS shapes and presets in WASM, graded domains, worker recovery/cancellation and STL/VTU export generation.
