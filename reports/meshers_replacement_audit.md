# Meshers replacement audit

The current experiment does not meet the replacement goal. Speed gains coexist with geometric-fidelity losses, inferior angle/size distributions in some cases, incomplete integration and unverified FEA accuracy. This audit inventories the known gaps across the current integration and marks untested paths explicitly. It does not claim every possible geometry and parameter combination has been benchmarked.

## Acceptance contract

Primary comparisons must use the same physical geometry and final triangle budget, or the same tetrahedron/DOF budget for solids. Counts and reference uncertainty must be reported, and any tolerance must be explicit. Surface quality requires angle and size distributions AND geometric fidelity. Simulation quality requires actual fedoo convergence, reactions, energy and solver checks. Fast/printing surfaces require fidelity, wall thickness and topology, even when element shape is unimportant. A failed baseline is a failure record, not a speed or accuracy win.

## Implementation order

### 1. Fix the comparison contract

Validate baseline topology and patch tests; implement explicit count budgets; retain failures and all costs. Add error-versus-count and time-to-error curves. Repair the MMG affine-patch baseline before using it for FEA comparisons.

### 2. Improve fast geometry

Test bounded one/two-step edge-root corrections and error-triggered refinement. Measure facet interiors and reverse distances. If the six-tetrahedron decomposition wastes the budget, prototype direct hexahedral extraction with consistent ambiguity handling and periodic caps.

### 3. Improve surface shape at fixed count

Implement periodic-orbit paired flips and tangential relocation, then count-neutral split/collapse redistribution. Optimize size distribution and angle tails together while rejecting geometric-error regressions. Start with the unit gyroid and graded Split-P tail defect.

### 4. Validate actual mechanics

For solids: affine patch, manufactured solution, compression, then six periodic homogenization loads on a refinement series. For shells: first generate genuine midsurfaces, then bending/membrane cases with the same thickness, formulation and BCs. Use a converged reference and solve-cost accounting.

### 5. Improve scaling and fill integration gaps

Use work-based worker selection, native field evaluation, deterministic block extraction, rolling caches and memory budgets. Then close anisotropy, chart seams/poles, sweeps, infill, density fitting, part types and generic-shape coverage. Release only the paths that meet their acceptance criteria.

## Complete gap registry

### 01 · Comparable budgets · P0 · Unproven

No exact output-triangle budget control. Prior comparisons have different counts; a smaller mesh is not automatically more efficient.

Evidence: surface_tradeoffs.json: MMGS relaxed versus quality meshers counts differ by about 0.8%, 23.9%, 2.9%, 38.2% for the four cases.

Improve: Add target_triangles and tolerance, then count-neutral flips/vertex motion and paired split-collapse redistribution. Choose reachable common budgets without changing the geometry.

Acceptance: Same final count for the primary comparison; report any unavoidable count tolerance explicitly and do not call it an exact-count win. Also compare time to equal FEA error and equal geometry error.

### 02 · Periodic gyroid surfaces · P0 · Measured loss

Quality mode loses to relaxed MMGS in both angle distribution and size uniformity.

Evidence: Unit gyroid: minimum angle 7.97 vs 15.53 degrees, p01 21.90 vs 34.24; area CV 0.460 vs 0.261. Counts 7572 vs 7508.

Improve: Current periodic polish moves vertices but freezes connectivity. Add orbit-paired edge flips first; then tangential relocation with a size term and strict geometry checks.

Acceptance: At the same count, improve angle tails and size distribution while keeping cap node/triangle matches, sampled geometry bounds and a time advantage.

### 03 · Graded gyroid surfaces · P0 · Measured loss

Quality mode improves the worst angle but loses the first-percentile angle and size uniformity.

Evidence: Graded gyroid: p01 25.30 vs MMGS 35.34; area CV 0.430 vs 0.280. Minimum angle is better: 14.67 vs 1.22.

Improve: Use a combined angle and sizing objective, then local connectivity edits. Track tails separately instead of optimizing only the single worst triangle.

Acceptance: Same-count quality curves and FEA error; no claim of overall quality superiority from one metric.

### 04 · Geometric tails after optimization · P0 · Measured loss

Quality mode can worsen the largest sampled geometric deviation.

Evidence: replacement_fidelity.json: graded Split-P max sampled distance 0.04118, VTK 0.02254, relaxed MMGS 0.01658. Reference-refinement sampled max 0.00509.

Improve: Audit collapse/flip/polish candidates near caps and high curvature. Bound facet-interior error, preserve features and reject moves exceeding the geometry budget.

Acceptance: Both high quantiles and worst sampled error improve at fixed count; add adaptive error sampling, feature coverage and self-intersection checks.

### 05 · Fast surface fidelity · P0 · Measured loss

Linear fast extraction loses geometric fidelity in all four tested cases.

Evidence: Bidirectional p95 errors are about 2.0-2.5 times VTK, despite similar counts. See replacement_fidelity.json.

Improve: Test one and two safeguarded edge-root correction steps; reject or refine cells with excessive facet error. Then test direct hexahedral extraction to reduce tetrahedral-diagonal overhead.

Acceptance: Same-count fast output must be faster and no worse in sampled distance tails, wall thickness, topology and volume error. Faster extraction alone fails the gate.

### 06 · Accurate surface fidelity · P0 · Measured loss

Exact edge roots do not ensure accurate facets at a fixed triangle count.

Evidence: Bidirectional p95 is worse than VTK for unit gyroid, graded gyroid and graded Split-P; graded Split-P max is 0.04098 vs VTK 0.02254.

Improve: Redistribute triangles by curvature and facet error; optimize projection and connectivity together. Avoid spending vertices on the six-tetrahedron grid decomposition.

Acceptance: Verify interiors and the reverse reference-to-output distance, not just vertices on the implicit zero set.

### 07 · Accurate surface speed · P1 · Measured loss

Exact-root mode is slower than raw VTK on the graded cases.

Evidence: Graded gyroid 0.100 vs 0.085 s; graded Split-P 0.235 vs 0.149 s. Counts are not identical.

Improve: Cache field/gradient evaluations, use safeguarded secant steps and analytic gradients, and retain early rejection. Profile extraction separately from metric computation and Python conversion.

Acceptance: Faster at a common output count AND no worse geometric fidelity.

### 08 · Raw volume generation · P1 · Measured loss

Zero-optimization tetrahedra are still slower than raw VTK mixed cells.

Evidence: fast_volume_data.json; VTK produces tetrahedra, hexahedra, wedges and pyramids, so output contracts differ.

Improve: Decide the required output contract first. For tetrahedra, optimize direct clipping and allocation; investigate a separate mixed-cell output only if microgen needs it.

Acceptance: Compare equivalent output types and validity. Keep the faster VTK path until meshers meets the actual user requirement.

### 09 · Volume quality recovery · P0 · Measured cost

A coarse Split-P quality failure can trigger global refinement and greatly inflate output.

Evidence: MESHERS_EXPERIMENT.md: 20,673 to 175,836 tetrahedra in the quality-recovery example.

Improve: Use local sliver repair, local refinement and constraint-aware relocation; expose recovery counts and timings. Avoid resolving one cap sliver by refining the entire background.

Acceptance: Meet quality and geometry thresholds at comparable tet/DOF budgets, then pass fedoo accuracy checks.

### 10 · Anisotropic cells · P1 · Measured failure

Tested anisotropic/repeated configurations fail quality or time out.

Evidence: experiments/tpms_meshers_results.json: anisotropic quality about 0.065-0.068 and a 90 s timeout.

Improve: Use physical-space edge metrics and balanced background spacing; separate chart distortion from element optimization.

Acceptance: Sweep cell aspect ratios and repeats, with signed volumes, periodic constraints and fedoo convergence.

### 11 · Coarse density grading · P1 · Measured failure

Some graded Split-P volumes fail the existing geometry gate.

Evidence: optimizer_experiments.json / .md: coarse graded Split-P error about 0.02163 versus tolerance 0.01.

Improve: Refine where field curvature and thickness require it, with a local error estimator; test gradient-aware sizing across the complete graded domain.

Acceptance: Pass geometry and solution error at bounded cost without relying on periodic tiling.

### 12 · Graded infill · P1 · Measured failure

A graded infill failed geometry acceptance and a finer callback run timed out.

Evidence: MESHERS_EXPERIMENT.md and experiments/tpms_meshers_results.json.

Improve: Compile or cache envelope-distance fields, add native interpolation, and use local refinement near envelope intersections.

Acceptance: Representative imported envelopes pass geometry/quality and run faster than the existing workflow.

### 13 · Density targets and sampled offsets · P1 · Integration gap

Direct surface presets reject density fitting and sampled offsets. Volume density fitting repeatedly meshes while searching.

Evidence: Tpms.generate_meshers_surface guards; generate_meshers density root search.

Improve: Calibrate volume fraction on reusable sampled fields, then validate on the final mesh; add native interpolation and a surface density-fitting adapter.

Acceptance: Preserve density semantics, quantify achieved density and total calibration cost, and restore inputs on failure.

### 14 · Full cylindrical wrap · P1 · Integration gap

Direct mapped output exists, but the wrap probe retained coincident unjoined seam points.

Evidence: curved_core_probe_data.json; full-wrap volume and surface returned meshes, so this is not proof of native inability.

Improve: Weld chart-seam equivalence classes with orientation-aware topology; carry rigid periodic transforms.

Acceptance: Closed conforming seam, no duplicate seam DOFs, and rotation-aware fedoo constraint checks.

### 15 · Spherical poles and collapsed axes · P1 · Measured failure

Full-sphere mapped volume probe failed with inverted background tetrahedra; a direct surface probe returned a mesh.

Evidence: curved_core_probe_data.json; regular sectors succeeded.

Improve: Use nonsingular charts or a Cartesian implicit representation near poles, then join interfaces consistently.

Acceptance: Positive volume elements, closure and no duplicate pole/seam nodes; validate solution convergence.

### 16 · Sweeps · P1 · Integration gap

No validated volume mapping is wired in; surface probe returned a mesh.

Evidence: curved_core_probe_data.json and Tpms._needs_parametric_clip.

Improve: Build a sweep chart with Jacobian checks, curvature limits and robust chart joining, or mesh a native Cartesian swept field.

Acceptance: Test curved paths, varying radius, closed loops and near self-contact with geometry and fedoo checks.

### 17 · Surface geometry coverage · P1 · Integration gap

The new surface method only accepts plain Cartesian Tpms sheets; curved types and Infill are rejected.

Evidence: Tpms.generate_meshers_surface type guard. Native curved surface probes demonstrate broader possibilities.

Improve: Connect validated implicit/chart surface paths to the adapter rather than labeling all these shapes unsupported by meshers.

Acceptance: Match the legacy API coverage and placement/labels for every enabled geometry.

### 18 · Surface part types · P0 · Integration gap

Skeletal surfaces and open zero-thickness midsurfaces are not exposed by the new band-sheet surface method.

Evidence: The new method fixes band=(-1,1) and has no type_part parameter.

Improve: Add one-sided inequalities and a single-isovalue/open-surface extractor; preserve boundary curves and normals.

Acceptance: This is required before meaningful shell FEA of a TPMS midsurface; do not substitute a closed solid boundary.

### 19 · Rectangular surface sampling · P1 · Integration gap

The surface adapter rejects unequal total grid counts; the native surface API accepts one cells value.

Evidence: triangles.rs extract_with_edge_refinement and Tpms.generate_meshers_surface.

Improve: Support per-axis cell counts and physical spacing in indexing, caching and cap extraction.

Acceptance: Anisotropic cells and repeats use the requested physical resolution without needless refinement.

### 20 · Automatic CPU selection · P0 · Measured regression

All-CPU defaults make small meshes slower.

Evidence: microgen_parallel_results.json: 0.890 s with one worker versus 1.541 s with 24. Earlier graded single mesh: 4.472 / 3.802 / 5.823 s at 1 / 2 / 4 workers.

Improve: Use a measured work-size threshold, persistent pools and parallel batches for small jobs. Do not equate maximum CPU use with minimum latency.

Acceptance: No small-case regression; demonstrate speedup and efficiency on larger single domains across worker counts.

### 21 · Python callback concurrency · P1 · Measured limit

Callbacks serialize through Python even though the native mesher releases the GIL.

Evidence: optimizer_experiments.md: four callback jobs take 2.219 / 2.204 / 2.286 s at 1 / 2 / 4 job threads.

Improve: Extend field compilation, native sampled-field interpolation and bulk evaluation. Consider process-based isolation only after measuring transfer and memory costs.

Acceptance: Benchmark both compiled and unavoidable callback fields; report effective worker use.

### 22 · Surface parallelism · P1 · Implementation gap

The direct surface extractor has serial grid sampling and extraction loops and no worker option.

Evidence: triangles.rs lines around extract_with_edge_refinement; Python generate_surface signature.

Improve: Parallel sample blocks and count output first; prefix sums allocate deterministic output ranges. Use local caches and deterministic seam reconciliation; color independent optimization patches.

Acceptance: Scaling on full graded domains, deterministic topology and periodicity, and no memory explosion or small-case slowdown.

### 23 · Grid size and memory · P1 · Scalability limit

Cells per axis are capped at 128. Extraction stores the full 3D point/value grid plus edge maps; large-job peak memory is not established.

Evidence: triangles.rs allocates n^3 points and values; native API validation; earlier scaling report lacks controlled peak RSS.

Improve: Use slab sampling/rolling edge caches or block extraction, reuse allocations, and add explicit memory budgets and failure cleanup.

Acceptance: Measure peak process memory and time versus physical volume and triangle count beyond unit cells; test graded, untileable geometry.

### 24 · Periodic tiling scope · P1 · Measured limit

Tiling cannot accelerate genuinely nonrepeating grading; mapped/intersection constraint metadata limits current tiling.

Evidence: TPMS_SCALING.md and periodic_tile.py.

Improve: Tile only verified identical directions; retain transforms and labels; improve direct whole-domain generation independently.

Acceptance: Whole-domain graded results must meet targets without presenting tiled periodic speedups as evidence.

### 25 · WASM latency · P2 · Measured gap

WASM quality generation remains slower than native; browser timing includes diagnostics and serialization.

Evidence: optimization_performance.md: gyroid 1.371 s WASM vs 0.649 s native in the recorded adapter experiment.

Improve: Profile phases, reuse buffers and compiled fields, use typed-array transfers and test supported SIMD/thread builds separately.

Acceptance: Record cold/warm module load, meshing, transfer and rendering independently on supported browsers.

### 26 · FEA accuracy at an element budget · P0 · Unproven

Good tetrahedron shape does not establish a sufficiently accurate structural response.

Evidence: New fedoo gyroid compression changes from stiffness 0.15916 at 7,213 tets to 0.12564 at 139,935; the last two levels still differ by about 3.3%.

Improve: Drive refinement by solution and geometry error, test tet10 or improved linear-tet sizing, and measure time to accuracy including solve cost.

Acceptance: Converged displacement, energy, reactions and effective stiffness; compare same tet/DOF budgets and total pipeline time against a valid baseline.

### 27 · Shell and periodic FEA coverage · P0 · Unverified

Current surface comparisons have no shell-FEA confirmation, and the new fedoo run is clamped compression plus manufactured solutions, not periodic homogenization.

Evidence: replacement_fedoo.json; older fedoo_homogenization results use a different specimen/platform.

Improve: Run genuine midsurface shell bending/membrane tests; run six bulk periodic strains, graded compression, thin features and stress convergence.

Acceptance: Preserve material, thickness, loading, formulation and BCs; independently audit energy and force balance. Do not infer shell validity from solid tetrahedra.

### 28 · MMG reference validity · P0 · Baseline defect

A failed baseline cannot establish comparative superiority.

Evidence: New 16-point MMG output has four same-side adjacent-tet face pairs; fedoo and independent affine tests both give L2 error about 5e-6. The 12-point MMG attempt timed out after 180 s.

Improve: Check original VTK triangulation, MMG output and boundary snapping separately. Find a conforming reference without dropping failed runs from the record.

Acceptance: Both baselines pass signed-volume, conforming-adjacency and affine-patch checks before comparison.

### 29 · Geometry guarantees · P0 · Validation gap

A p95 implicit residual is neither a distance bound nor a topology/thickness guarantee.

Evidence: New bidirectional distance checks reveal graded Split-P tail regression missed by average metrics.

Improve: Use converged independent references, bidirectional samples with adaptive refinement, topology, intersections, solid volume and minimum wall-thickness checks.

Acceptance: State sampling limits. Require accuracy margins larger than reference uncertainty, and certify tails before calling a mesh suitable for printing.

### 30 · Existing mesh remeshing · P1 · Integration gap

MMG also adapts imported meshes and preserves boundary entities; an implicit generator does not replace that workflow.

Evidence: microgen.remesh.remesh_keeping_boundaries_for_fem and external.Mmg APIs.

Improve: Inventory boundary groups, size fields, existing-mesh adaptation and metadata contracts; implement equivalents or retain the dependency for those paths.

Acceptance: Every public remeshing operation has tested equivalent behavior before removing MMG.

### 31 · Other microgen geometries · P1 · Unverified coverage

Primitive shapes, strut lattices, spinodoids, booleans and imported polyhedra have not been shown to improve through the new TPMS adapter.

Evidence: microgen.shape public class inventory; Shape.generate_surface_mesh remains its own workflow.

Improve: Add a common implicit-generation adapter and representative per-family benchmarks; preserve feature edges, disconnected components and material labels.

Acceptance: Classify each path as improved, unchanged, failing or untested. CAD/B-rep operations remain separate capabilities.

### 32 · Experimental APIs and presets · P1 · Release gap

Direct surface support requires an unreleased feature build; legacy surface API still runs VTK, and preset names are not acceptance guarantees.

Evidence: MESHERS_EXPERIMENT.md and the new generate_meshers_surface method.

Improve: Stabilize a quality/error/budget contract, ship supported wheels, integrate normal API routing only after acceptance tests, and display achieved diagnostics.

Acceptance: Cross-platform packaging, cancellation/failure behavior, metadata and reproducible versioned benchmarks. No silent quality downgrade or hidden backend choice.

## Fedoo findings

Five meshers gyroid refinement levels pass affine patch and manufactured-solution checks. Compression stiffness is not yet converged: the final 28-to-36 sample refinement changes stiffness by about 3.3%. The 16-sample stiffness is about 15% above the finest tested value, which is itself not an exact reference. This is clamped compression with free side surfaces, not periodic homogenization. No shell-FEA result is claimed.

The 16-sample MMG baseline failed both fedoo and independent affine tests; four internal faces have same-side adjacent tetrahedra. Double-precision reruns did not fix it. The 20-sample fedoo patch also failed; its topology has not yet been diagnosed. The 12-sample attempt timed out after 180 seconds. These outputs do not provide a valid comparative FEA reference.