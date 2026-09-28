"""Render the table-driven extraction experiment inside the limitations report."""

import html
import json
import statistics
from pathlib import Path


def section():
    data = json.loads((Path(__file__).parent / "table_surface_data.json").read_text())
    cases = list(data["selected_resolutions"])
    rows, quality_rows, chart = [], [], []
    timings = {}
    modes = ("vtk", "meshers_fast", "meshers_linear")
    for index, case in enumerate(cases):
        groups = [
            [r for r in data["runs"] if r["case"] == case and r["mode"] == mode]
            for mode in modes
        ]
        assert all(
            len(group) == 3 and all(r["status"] == "ok" for r in group)
            for group in groups
        )
        for group in groups:
            assert all(r["open_edges"] == 0 for r in group)
            assert all(
                not any(r["periodic_node_mismatch"] + r["periodic_triangle_mismatch"])
                for r in group
            )
        times = [
            statistics.median(r["generation_seconds"] for r in group)
            for group in groups
        ]
        timings[case] = times
        counts = [group[0]["triangle_count"] for group in groups]
        gap = counts[2] / counts[0] - 1
        name = html.escape(case.replace("_", " ").replace("unit", "1 cell"))
        rows.append(
            f"<tr><th>{name}</th><td>{times[0]:.3f} s</td><td>{times[1]:.3f} s</td><td>{times[2]:.3f} s</td><td>{times[0] / times[2]:.2f}×</td><td>{counts[0]:,} / {counts[2]:,}</td><td>{gap:+.1%}</td></tr>"
        )
        for mode, group in zip(modes, groups):
            r = group[0]
            quality_rows.append(
                f"<tr><th>{name}</th><td>{mode}</td><td>{r['angle_p01_degrees']:.2f}°</td><td>{r['area_cv']:.3f}</td><td>{r['equivalent_diameter_p01']:.4f} / {r['equivalent_diameter_median']:.4f} / {r['equivalent_diameter_p99']:.4f}</td><td>{r['implicit_vertex_residual_p95']:.3g}</td></tr>"
            )
        y = 40 + index * 44
        ratio = times[2] / times[0]
        chart.append(
            f'<text x="10" y="{y + 15}" class="row-label">{name}</text><rect x="225" y="{y}" width="{ratio * 500:.1f}" height="24" rx="4" fill="var(--blue)"/><text x="{235 + ratio * 500:.1f}" y="{y + 17}" class="ratio">{ratio:.2f}×</text>'
        )
    wins = sum(t[2] < t[0] for t in timings.values())
    speedups = [t[0] / t[2] for t in timings.values()]
    scale = []
    for size in (2, 3, 5, 7):
        t = timings[f"graded_gyroid_{size}"]
        scale.append(
            f"<tr><th>{size}³ = {size**3} cells</th><td>{t[0]:.3f} s</td><td>{t[2]:.3f} s</td><td>{1000 * t[0] / size**3:.2f}</td><td>{1000 * t[2] / size**3:.2f}</td></tr>"
        )
    return f'''<section id="table-extractor">
<div class="section-tag">New experiment · {html.escape(data["measured_at_utc"][:10])}</div>
<h2>Table-driven extraction beats raw VTK in {wins} of {len(cases)} cases</h2>
<p>The Rust linear prototype now uses a 16-case tetrahedron table. It avoids implicit-wall polygon sorting and analytic normals, checks box boundaries directly, and avoids adjacency setup for zero-pass quality calls. It still uses six tetrahedra per voxel. This is not a port of VTK's hexahedron clipping tables or Flying Edges.</p>
<div class="note">Measured speedup over microgen's complete raw VTK surface path: {min(speedups):.2f}–{max(speedups):.2f}× at similar triangle counts. This proves that meshers can beat that workflow on these cases. It does not establish equal geometric accuracy or performance against VTK Flying Edges.</div>
<div class="table-wrap"><table><thead><tr><th>Case</th><th>VTK</th><th>Meshers exact, no polish</th><th>Meshers table linear</th><th>Linear speedup</th><th>Triangles VTK / linear</th><th>Count gap</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>
<div class="chart-card"><svg class="plot" role="img" aria-label="Linear meshers generation time divided by VTK time, lower is faster" viewBox="0 0 850 {60 + 44 * len(cases)}"><text x="225" y="22" class="axis-title">Meshers linear / microgen VTK time, lower is faster</text><line x1="725" y1="30" x2="725" y2="{40 + 44 * len(cases)}" stroke="var(--vtk)" stroke-dasharray="5 5"/><text x="730" y="22" class="tick">VTK = 1×</text>{"".join(chart)}</svg></div>
<h3>Scaling with density grading</h3>
<p>These domains are meshed in full. No unit cell is replicated. Grading uses the same spatial slope, so the thickness range grows with the domain. The 2³ case uses 12 meshers points per cell; 3³, 5³ and 7³ use 11. VTK uses 16 throughout. Timings include field construction and sampling, extraction, Python output and meshers' built-in checks. Imports, process startup, and the additional quality and topology measurements are excluded. Memory scaling was not measured.</p>
<div class="table-wrap"><table><thead><tr><th>Domain</th><th>VTK</th><th>Table linear</th><th>VTK ms / cell</th><th>Linear ms / cell</th></tr></thead><tbody>{"".join(scale)}</tbody></table></div>
<h3>What the faster path gives up</h3>
<p>Linear intersections retain the earlier prototype's geometric error. Triangle counts, angles and size distribution alone do not establish printable fidelity. Exact intersections remain available and the default quality settings are unchanged. Exact roots constrain vertices to the implicit surface, but neither vertex residual nor closure measures facet-interior error, minimum wall thickness or self-intersection.</p>
<p>Every timed output had zero open edges and matching nodes and triangles on the requested periodic faces. Fully graded domains have no periodic axes. Metrics below include caps. Angle is the 1st percentile of minimum triangle angle; area CV is standard deviation divided by mean. Equivalent-area diameter describes triangle size in model units. Residual is the 95th percentile of ||f|−1| at interior surface vertices, a dimensionless value rather than a physical distance.</p>
<details><summary>Quality, triangle size distribution and geometric residual</summary><div class="table-wrap"><table><thead><tr><th>Case</th><th>Path</th><th>Angle p01</th><th>Area CV</th><th>Diameter p01 / median / p99</th><th>Vertex residual p95</th></tr></thead><tbody>{"".join(quality_rows)}</tbody></table></div></details>
<p>Each path ran three times in fresh processes, with order alternated. The first six cases reuse the previously selected count-matched resolutions. Larger graded cases retain the 3³ sampling density; their actual count differences are reported above. These results come from one host and do not establish a universal speed guarantee.</p>
<p>Microgen uses two <code>vtkTableBasedClipDataSet</code> calls, then surface extraction, cleaning and triangulation. <a href="https://vtk.org/doc/nightly/html/classvtkTableBasedClipDataSet.html">VTK's clipping documentation</a> describes its case-table method. <a href="https://vtk.org/doc/nightly/html/classvtkFlyingEdges3D.html">Flying Edges</a> provides a separate model for structured-grid surface extraction with single edge intersections and preallocated output. A direct hexahedral surface extractor is a further experiment, especially for reducing geometric error at a fixed triangle budget.</p>
<p class="muted">Meshers source revision <code>{html.escape(data['meshers_commit'])}</code>. Reproduce with <a href="table_surface_benchmark.py">the benchmark runner</a>, using the experimental meshers package and examples on PYTHONPATH. <a href="table_surface_data.json">All 72 measurements and quality metrics</a>. The earlier results below are retained as historical baselines.</p>
</section>'''
