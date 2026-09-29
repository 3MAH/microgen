"""Known replacement gaps, with evidence and concrete next experiments."""

import json
from pathlib import Path

root = Path(__file__).resolve().parent
rows = []


def add(id, status, area, problem, evidence, improvement, acceptance, priority="P1"):
    rows.append(
        dict(
            id=id,
            status=status,
            area=area,
            problem=problem,
            evidence=evidence,
            improvement=improvement,
            acceptance=acceptance,
            priority=priority,
        )
    )


add(
    "01",
    "Unproven",
    "Comparable budgets",
    "No exact output-triangle budget control. Prior comparisons have different counts; a smaller mesh is not automatically more efficient.",
    "surface_tradeoffs.json: MMGS relaxed versus quality meshers counts differ by about 0.8%, 23.9%, 2.9%, 38.2% for the four cases.",
    "Add target_triangles and tolerance, then count-neutral flips/vertex motion and paired split-collapse redistribution. Choose reachable common budgets without changing the geometry.",
    "Same final count for the primary comparison; report any unavoidable count tolerance explicitly and do not call it an exact-count win. Also compare time to equal FEA error and equal geometry error.",
    "P0",
)
add(
    "02",
    "Measured loss",
    "Periodic gyroid surfaces",
    "Quality mode loses to relaxed MMGS in both angle distribution and size uniformity.",
    "Unit gyroid: minimum angle 7.97 vs 15.53 degrees, p01 21.90 vs 34.24; area CV 0.460 vs 0.261. Counts 7572 vs 7508.",
    "Current periodic polish moves vertices but freezes connectivity. Add orbit-paired edge flips first; then tangential relocation with a size term and strict geometry checks.",
    "At the same count, improve angle tails and size distribution while keeping cap node/triangle matches, sampled geometry bounds and a time advantage.",
    "P0",
)
add(
    "03",
    "Measured loss",
    "Graded gyroid surfaces",
    "Quality mode improves the worst angle but loses the first-percentile angle and size uniformity.",
    "Graded gyroid: p01 25.30 vs MMGS 35.34; area CV 0.430 vs 0.280. Minimum angle is better: 14.67 vs 1.22.",
    "Use a combined angle and sizing objective, then local connectivity edits. Track tails separately instead of optimizing only the single worst triangle.",
    "Same-count quality curves and FEA error; no claim of overall quality superiority from one metric.",
    "P0",
)
add(
    "04",
    "Measured loss",
    "Geometric tails after optimization",
    "Quality mode can worsen the largest sampled geometric deviation.",
    "replacement_fidelity.json: graded Split-P max sampled distance 0.04118, VTK 0.02254, relaxed MMGS 0.01658. Reference-refinement sampled max 0.00509.",
    "Audit collapse/flip/polish candidates near caps and high curvature. Bound facet-interior error, preserve features and reject moves exceeding the geometry budget.",
    "Both high quantiles and worst sampled error improve at fixed count; add adaptive error sampling, feature coverage and self-intersection checks.",
    "P0",
)
add(
    "05",
    "Measured loss",
    "Fast surface fidelity",
    "Linear fast extraction loses geometric fidelity in all four tested cases.",
    "Bidirectional p95 errors are about 2.0-2.5 times VTK, despite similar counts. See replacement_fidelity.json.",
    "Test one and two safeguarded edge-root correction steps; reject or refine cells with excessive facet error. Then test direct hexahedral extraction to reduce tetrahedral-diagonal overhead.",
    "Same-count fast output must be faster and no worse in sampled distance tails, wall thickness, topology and volume error. Faster extraction alone fails the gate.",
    "P0",
)
add(
    "06",
    "Measured loss",
    "Accurate surface fidelity",
    "Exact edge roots do not ensure accurate facets at a fixed triangle count.",
    "Bidirectional p95 is worse than VTK for unit gyroid, graded gyroid and graded Split-P; graded Split-P max is 0.04098 vs VTK 0.02254.",
    "Redistribute triangles by curvature and facet error; optimize projection and connectivity together. Avoid spending vertices on the six-tetrahedron grid decomposition.",
    "Verify interiors and the reverse reference-to-output distance, not just vertices on the implicit zero set.",
    "P0",
)
add(
    "07",
    "Measured loss",
    "Accurate surface speed",
    "Exact-root mode is slower than raw VTK on the graded cases.",
    "Graded gyroid 0.100 vs 0.085 s; graded Split-P 0.235 vs 0.149 s. Counts are not identical.",
    "Cache field/gradient evaluations, use safeguarded secant steps and analytic gradients, and retain early rejection. Profile extraction separately from metric computation and Python conversion.",
    "Faster at a common output count AND no worse geometric fidelity.",
    "P1",
)
add(
    "08",
    "Measured loss",
    "Raw volume generation",
    "Zero-optimization tetrahedra are still slower than raw VTK mixed cells.",
    "fast_volume_data.json; VTK produces tetrahedra, hexahedra, wedges and pyramids, so output contracts differ.",
    "Decide the required output contract first. For tetrahedra, optimize direct clipping and allocation; investigate a separate mixed-cell output only if microgen needs it.",
    "Compare equivalent output types and validity. Keep the faster VTK path until meshers meets the actual user requirement.",
    "P1",
)
add(
    "09",
    "Measured cost",
    "Volume quality recovery",
    "A coarse Split-P quality failure can trigger global refinement and greatly inflate output.",
    "MESHERS_EXPERIMENT.md: 20,673 to 175,836 tetrahedra in the quality-recovery example.",
    "Use local sliver repair, local refinement and constraint-aware relocation; expose recovery counts and timings. Avoid resolving one cap sliver by refining the entire background.",
    "Meet quality and geometry thresholds at comparable tet/DOF budgets, then pass fedoo accuracy checks.",
    "P0",
)
add(
    "10",
    "Measured failure",
    "Anisotropic cells",
    "Tested anisotropic/repeated configurations fail quality or time out.",
    "experiments/tpms_meshers_results.json: anisotropic quality about 0.065-0.068 and a 90 s timeout.",
    "Use physical-space edge metrics and balanced background spacing; separate chart distortion from element optimization.",
    "Sweep cell aspect ratios and repeats, with signed volumes, periodic constraints and fedoo convergence.",
    "P1",
)
add(
    "11",
    "Measured failure",
    "Coarse density grading",
    "Some graded Split-P volumes fail the existing geometry gate.",
    "optimizer_experiments.json / .md: coarse graded Split-P error about 0.02163 versus tolerance 0.01.",
    "Refine where field curvature and thickness require it, with a local error estimator; test gradient-aware sizing across the complete graded domain.",
    "Pass geometry and solution error at bounded cost without relying on periodic tiling.",
    "P1",
)
add(
    "12",
    "Measured failure",
    "Graded infill",
    "A graded infill failed geometry acceptance and a finer callback run timed out.",
    "MESHERS_EXPERIMENT.md and experiments/tpms_meshers_results.json.",
    "Compile or cache envelope-distance fields, add native interpolation, and use local refinement near envelope intersections.",
    "Representative imported envelopes pass geometry/quality and run faster than the existing workflow.",
    "P1",
)
add(
    "13",
    "Integration gap",
    "Density targets and sampled offsets",
    "Direct surface presets reject density fitting and sampled offsets. Volume density fitting repeatedly meshes while searching.",
    "Tpms.generate_meshers_surface guards; generate_meshers density root search.",
    "Calibrate volume fraction on reusable sampled fields, then validate on the final mesh; add native interpolation and a surface density-fitting adapter.",
    "Preserve density semantics, quantify achieved density and total calibration cost, and restore inputs on failure.",
    "P1",
)
add(
    "14",
    "Integration gap",
    "Full cylindrical wrap",
    "Direct mapped output exists, but the wrap probe retained coincident unjoined seam points.",
    "curved_core_probe_data.json; full-wrap volume and surface returned meshes, so this is not proof of native inability.",
    "Weld chart-seam equivalence classes with orientation-aware topology; carry rigid periodic transforms.",
    "Closed conforming seam, no duplicate seam DOFs, and rotation-aware fedoo constraint checks.",
    "P1",
)
add(
    "15",
    "Measured failure",
    "Spherical poles and collapsed axes",
    "Full-sphere mapped volume probe failed with inverted background tetrahedra; a direct surface probe returned a mesh.",
    "curved_core_probe_data.json; regular sectors succeeded.",
    "Use nonsingular charts or a Cartesian implicit representation near poles, then join interfaces consistently.",
    "Positive volume elements, closure and no duplicate pole/seam nodes; validate solution convergence.",
    "P1",
)
add(
    "16",
    "Integration gap",
    "Sweeps",
    "No validated volume mapping is wired in; surface probe returned a mesh.",
    "curved_core_probe_data.json and Tpms._needs_parametric_clip.",
    "Build a sweep chart with Jacobian checks, curvature limits and robust chart joining, or mesh a native Cartesian swept field.",
    "Test curved paths, varying radius, closed loops and near self-contact with geometry and fedoo checks.",
    "P1",
)
add(
    "17",
    "Integration gap",
    "Surface geometry coverage",
    "The new surface method only accepts plain Cartesian Tpms sheets; curved types and Infill are rejected.",
    "Tpms.generate_meshers_surface type guard. Native curved surface probes demonstrate broader possibilities.",
    "Connect validated implicit/chart surface paths to the adapter rather than labeling all these shapes unsupported by meshers.",
    "Match the legacy API coverage and placement/labels for every enabled geometry.",
    "P1",
)
add(
    "18",
    "Integration gap",
    "Surface part types",
    "Skeletal surfaces and open zero-thickness midsurfaces are not exposed by the new band-sheet surface method.",
    "The new method fixes band=(-1,1) and has no type_part parameter.",
    "Add one-sided inequalities and a single-isovalue/open-surface extractor; preserve boundary curves and normals.",
    "This is required before meaningful shell FEA of a TPMS midsurface; do not substitute a closed solid boundary.",
    "P0",
)
add(
    "19",
    "Integration gap",
    "Rectangular surface sampling",
    "The surface adapter rejects unequal total grid counts; the native surface API accepts one cells value.",
    "triangles.rs extract_with_edge_refinement and Tpms.generate_meshers_surface.",
    "Support per-axis cell counts and physical spacing in indexing, caching and cap extraction.",
    "Anisotropic cells and repeats use the requested physical resolution without needless refinement.",
    "P1",
)
add(
    "20",
    "Measured regression",
    "Automatic CPU selection",
    "All-CPU defaults make small meshes slower.",
    "microgen_parallel_results.json: 0.890 s with one worker versus 1.541 s with 24. Earlier graded single mesh: 4.472 / 3.802 / 5.823 s at 1 / 2 / 4 workers.",
    "Use a measured work-size threshold, persistent pools and parallel batches for small jobs. Do not equate maximum CPU use with minimum latency.",
    "No small-case regression; demonstrate speedup and efficiency on larger single domains across worker counts.",
    "P0",
)
add(
    "21",
    "Measured limit",
    "Python callback concurrency",
    "Callbacks serialize through Python even though the native mesher releases the GIL.",
    "optimizer_experiments.md: four callback jobs take 2.219 / 2.204 / 2.286 s at 1 / 2 / 4 job threads.",
    "Extend field compilation, native sampled-field interpolation and bulk evaluation. Consider process-based isolation only after measuring transfer and memory costs.",
    "Benchmark both compiled and unavoidable callback fields; report effective worker use.",
    "P1",
)
add(
    "22",
    "Implementation gap",
    "Surface parallelism",
    "The direct surface extractor has serial grid sampling and extraction loops and no worker option.",
    "triangles.rs lines around extract_with_edge_refinement; Python generate_surface signature.",
    "Parallel sample blocks and count output first; prefix sums allocate deterministic output ranges. Use local caches and deterministic seam reconciliation; color independent optimization patches.",
    "Scaling on full graded domains, deterministic topology and periodicity, and no memory explosion or small-case slowdown.",
    "P1",
)
add(
    "23",
    "Scalability limit",
    "Grid size and memory",
    "Cells per axis are capped at 128. Extraction stores the full 3D point/value grid plus edge maps; large-job peak memory is not established.",
    "triangles.rs allocates n^3 points and values; native API validation; earlier scaling report lacks controlled peak RSS.",
    "Use slab sampling/rolling edge caches or block extraction, reuse allocations, and add explicit memory budgets and failure cleanup.",
    "Measure peak process memory and time versus physical volume and triangle count beyond unit cells; test graded, untileable geometry.",
    "P1",
)
add(
    "24",
    "Measured limit",
    "Periodic tiling scope",
    "Tiling cannot accelerate genuinely nonrepeating grading; mapped/intersection constraint metadata limits current tiling.",
    "TPMS_SCALING.md and periodic_tile.py.",
    "Tile only verified identical directions; retain transforms and labels; improve direct whole-domain generation independently.",
    "Whole-domain graded results must meet targets without presenting tiled periodic speedups as evidence.",
    "P1",
)
add(
    "25",
    "Measured gap",
    "WASM latency",
    "WASM quality generation remains slower than native; browser timing includes diagnostics and serialization.",
    "optimization_performance.md: gyroid 1.371 s WASM vs 0.649 s native in the recorded adapter experiment.",
    "Profile phases, reuse buffers and compiled fields, use typed-array transfers and test supported SIMD/thread builds separately.",
    "Record cold/warm module load, meshing, transfer and rendering independently on supported browsers.",
    "P2",
)
add(
    "26",
    "Unproven",
    "FEA accuracy at an element budget",
    "Good tetrahedron shape does not establish a sufficiently accurate structural response.",
    "New fedoo gyroid compression changes from stiffness 0.15916 at 7,213 tets to 0.12564 at 139,935; the last two levels still differ by about 3.3%.",
    "Drive refinement by solution and geometry error, test tet10 or improved linear-tet sizing, and measure time to accuracy including solve cost.",
    "Converged displacement, energy, reactions and effective stiffness; compare same tet/DOF budgets and total pipeline time against a valid baseline.",
    "P0",
)
add(
    "27",
    "Unverified",
    "Shell and periodic FEA coverage",
    "Current surface comparisons have no shell-FEA confirmation, and the new fedoo run is clamped compression plus manufactured solutions, not periodic homogenization.",
    "replacement_fedoo.json; older fedoo_homogenization results use a different specimen/platform.",
    "Run genuine midsurface shell bending/membrane tests; run six bulk periodic strains, graded compression, thin features and stress convergence.",
    "Preserve material, thickness, loading, formulation and BCs; independently audit energy and force balance. Do not infer shell validity from solid tetrahedra.",
    "P0",
)
add(
    "28",
    "Baseline defect",
    "MMG reference validity",
    "A failed baseline cannot establish comparative superiority.",
    "New 16-point MMG output has four same-side adjacent-tet face pairs; fedoo and independent affine tests both give L2 error about 5e-6. The 12-point MMG attempt timed out after 180 s.",
    "Check original VTK triangulation, MMG output and boundary snapping separately. Find a conforming reference without dropping failed runs from the record.",
    "Both baselines pass signed-volume, conforming-adjacency and affine-patch checks before comparison.",
    "P0",
)
add(
    "29",
    "Validation gap",
    "Geometry guarantees",
    "A p95 implicit residual is neither a distance bound nor a topology/thickness guarantee.",
    "New bidirectional distance checks reveal graded Split-P tail regression missed by average metrics.",
    "Use converged independent references, bidirectional samples with adaptive refinement, topology, intersections, solid volume and minimum wall-thickness checks.",
    "State sampling limits. Require accuracy margins larger than reference uncertainty, and certify tails before calling a mesh suitable for printing.",
    "P0",
)
add(
    "30",
    "Integration gap",
    "Existing mesh remeshing",
    "MMG also adapts imported meshes and preserves boundary entities; an implicit generator does not replace that workflow.",
    "microgen.remesh.remesh_keeping_boundaries_for_fem and external.Mmg APIs.",
    "Inventory boundary groups, size fields, existing-mesh adaptation and metadata contracts; implement equivalents or retain the dependency for those paths.",
    "Every public remeshing operation has tested equivalent behavior before removing MMG.",
    "P1",
)
add(
    "31",
    "Unverified coverage",
    "Other microgen geometries",
    "Primitive shapes, strut lattices, spinodoids, booleans and imported polyhedra have not been shown to improve through the new TPMS adapter.",
    "microgen.shape public class inventory; Shape.generate_surface_mesh remains its own workflow.",
    "Add a common implicit-generation adapter and representative per-family benchmarks; preserve feature edges, disconnected components and material labels.",
    "Classify each path as improved, unchanged, failing or untested. CAD/B-rep operations remain separate capabilities.",
    "P1",
)
add(
    "32",
    "Release gap",
    "Experimental APIs and presets",
    "Direct surface support requires an unreleased feature build; legacy surface API still runs VTK, and preset names are not acceptance guarantees.",
    "MESHERS_EXPERIMENT.md and the new generate_meshers_surface method.",
    "Stabilize a quality/error/budget contract, ship supported wheels, integrate normal API routing only after acceptance tests, and display achieved diagnostics.",
    "Cross-platform packaging, cancellation/failure behavior, metadata and reproducible versioned benchmarks. No silent quality downgrade or hidden backend choice.",
    "P1",
)
(root / "replacement_gaps.json").write_text(json.dumps(rows, indent=2))
print(len(rows), "gaps recorded")
