"""Render recorded surface benchmarks without external browser dependencies."""

import json
import math
import statistics
from pathlib import Path

root = Path(__file__).resolve().parent
data = json.loads((root / "surface_tradeoffs.json").read_text())
modes = ["vtk", "mmgs", "mmgs_count", "fast", "accurate", "quality"]
names = {
    "vtk": "VTK raw",
    "mmgs": "VTK + MMGS tight",
    "mmgs_count": "VTK + MMGS relaxed",
    "fast": "Meshers fast",
    "accurate": "Meshers accurate",
    "quality": "Meshers quality",
}
colors = ["#9baabc", "#f08080", "#ffa94d", "#59d5d0", "#81a4ff", "#c49bff"]
cases = list(dict.fromkeys(r["case"] for r in data["runs"]))
parts = []
for case in cases:
    rows = []
    for mode in modes:
        records = [
            r
            for r in data["runs"]
            if r["case"] == case and r["mode"] == mode and r["status"] == "ok"
        ]
        if not records:
            continue
        row = dict(records[0])
        row["seconds"] = statistics.median(r["seconds"] for r in records)
        row["n"] = len(records)
        row["range"] = (
            f"{min(r['seconds'] for r in records):.3f}–{max(r['seconds'] for r in records):.3f}"
        )
        rows.append(row)
    svg = [
        '<svg viewBox="0 0 640 280" role="img" aria-label="Generation time versus first percentile triangle angle">'
    ]
    lo = min(math.log10(r["seconds"]) for r in rows) - 0.15
    hi = max(math.log10(r["seconds"]) for r in rows) + 0.15
    for q in (0, 10, 20, 30):
        y = 230 - q * 6
        svg.append(
            f'<path d="M60 {y}H610" stroke="#344352"/><text x="20" y="{y + 4}">{q}°</text>'
        )
    for t in (0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100):
        if lo <= math.log10(t) <= hi:
            x = 60 + (math.log10(t) - lo) / (hi - lo) * 550
            svg.append(f'<text x="{x}" y="253" text-anchor="middle">{t:g}s</text>')
    for r in rows:
        x = 60 + (math.log10(r["seconds"]) - lo) / (hi - lo) * 550
        y = 230 - r["angle_p01_degrees"] * 6
        color = colors[modes.index(r["mode"])]
        svg.append(
            f'<circle cx="{x}" cy="{y}" r="7" fill="{color}"><title>{names[r["mode"]]}: {r["seconds"]:.4f}s, {r["angle_p01_degrees"]:.2f}°</title></circle>'
        )
    svg.append(
        '<text x="320" y="277" text-anchor="middle">Time, logarithmic scale · faster ← | ↑ better angles</text></svg>'
    )
    table = []
    for r in rows:
        mismatches = r["periodic_node_mismatch"] + r["periodic_triangle_mismatch"]
        periodic = (
            "not requested"
            if not mismatches
            else ("matched" if not any(mismatches) else "FAILED")
        )
        if any(0 in pair for pair in r.get("periodic_cap_triangle_counts", [])):
            periodic = "MISSING CAP"
        gap = 100 * (r["triangle_count"] / r["target_triangles"] - 1)
        diameter = " / ".join(
            f"{r[k]:.4f}"
            for k in (
                "equivalent_diameter_p01",
                "equivalent_diameter_median",
                "equivalent_diameter_p99",
            )
        )
        table.append(
            f"<tr><th>{names[r['mode']]}</th><td>{r['seconds']:.3f}<small>{r['range']}</small></td><td>{r['triangle_count']:,}<small>{gap:+.1f}% vs VTK</small></td><td>{r['minimum_angle_degrees']:.2f} / {r['angle_p01_degrees']:.2f}</td><td>{r['area_cv']:.3f}</td><td>{diameter}</td><td>{r['sampled_distance_p95']:.5f} / {r['sampled_distance_max']:.5f}</td><td>{r['open_edges']} / {r['nonmanifold_edges']}</td><td>{periodic}</td></tr>"
        )
    parts.append(
        f'<section><h2>{case.replace("_", " ")}</h2>{"".join(svg)}<div class="scroll"><table><tr><th>Path</th><th>Seconds</th><th>Triangles</th><th>Angle min / p01, °</th><th>Area CV</th><th>Diameter p01 / median / p99</th><th>Distance proxy p95 / max</th><th>Open / nonmanifold edges</th><th>Periodic faces</th></tr>{"".join(table)}</table></div></section>'
    )
legend = " ".join(
    f'<span style="color:{c}">● {names[m]}</span>' for m, c in zip(modes, colors)
)
page = (
    """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>TPMS surface quality and performance</title><style>
body{background:#101822;color:#e4edf6;font:16px/1.55 system-ui;margin:0;padding:36px;max-width:1450px;margin:auto}h1{font-size:38px;line-height:1.15}h2{font-size:23px}p{max-width:1050px}section{margin:32px 0;padding:24px;background:#182431;border:1px solid #33475a;border-radius:14px}a{color:#7fd9ec}small{display:block;color:#a9bbce}table{border-collapse:collapse;font-size:13px;width:100%}td,th{padding:10px;text-align:left;border-bottom:1px solid #344352}th{color:#c7d7e7}.scroll{overflow-x:auto}svg{width:min(100%,700px);display:block;margin:15px 0}svg text{fill:#b8c9d9;font-size:12px}.legend span{display:inline-block;margin-right:18px}pre{padding:18px;background:#0e151e;overflow:auto}code{color:#8be0d5}.note{border-left:4px solid #f6be67;padding:12px 20px;background:#263040}</style>
<h1>TPMS surfaces: time, quality and fidelity</h1><p>Native Windows measurements · 29 September 2026 · unit-cell gyroid and Split-P, plus density-graded 2³ domains. Each path has three timed runs. Timings show medians and ranges. Surface quality is measured on all triangles, including clipping caps.</p>
<p class="note">These results compare specific MMGS settings, not the best possible MMGS configuration. MMGS starts from microgen's VTK surface and cannot recover the original implicit geometry exactly. Periodic matching is measured, not assumed. A quality preset is an optimization budget, not an FEA or printing certification.</p>
<h2>Choose the work budget</h2><p><b>fast</b> uses linear edge intersections and no triangle optimization. <b>accurate</b> solves implicit edge intersections but does not improve triangle shapes. <b>quality</b> adds paired polishing for periodic surfaces or topology improvements for nonperiodic surfaces. Exact vertices do not make curved triangle interiors exact.</p>
<pre><code>surface = tpms.generate_meshers_surface(
    optimization="quality", periodic=(True, True, True)
)
# Explicit controls override the preset.
surface = tpms.generate_meshers_surface(
    optimization="quality", periodic=(True, True, True), polish_passes=6
)</code></pre><p>This new API requires the experimental meshers build. It supports plain sheet TPMS with scalar or callable offsets and equal grid counts on all axes. Density fitting, curved charts and skeletal parts are not wired into this surface method. The legacy surface API remains available. Callable sheet offsets must remain positive throughout the domain.</p>
"""
    + f'<div class="legend">{legend}</div>'
    + "".join(parts)
    + """
<h2>How to read the comparison</h2><p>Triangle budgets are approximate. Meshers uses 12 samples per cell, except unit Split-P with 13. VTK uses 16. Quality optimization can change topology or retry the background grid. MMGS requested edge size is tuned separately toward the VTK count, with up to five trials for tight tolerance and eight for relaxed tolerance. Tight hausd is 0.001; relaxed hausd is 0.01, in model units. Tuning cost is excluded and all actual output counts are shown. MMGS can retain more triangles despite increasing the requested edge size.</p>
<p>Timing includes geometry setup and generation, plus file writing, executable launch and output loading for MMGS. It excludes module imports and metric evaluation. Meshers uses compiled native fields. Individual paths are run sequentially to avoid contention; the experimental surface extractor does not expose native worker selection. The raw JSON retains any preliminary native-only measurements separately.</p>
<p>Angles are the minimum angle in each triangle. Area CV is standard deviation divided by mean, with lower values indicating more uniform areas. Equivalent-area diameter is 2√(area/π). The distance proxy samples each wall triangle's centroid and three edge midpoints and evaluates | |f|−1 | / |∇f| against the original normalized implicit field. Caps are excluded. This is a first-order distance estimate, not a Hausdorff bound; near stationary points the estimate can be large. No minimum-wall-thickness, self-intersection or FEA convergence certification is implied.</p>
<p>Periodic validation compares both nodes and cap triangles on the original requested box planes, using a 1e-8 plane tolerance and rounding transverse coordinates to 1e-8 model units. Empty cap comparisons cannot pass. Grading along all three axes disables periodic constraints. Zero open and nonmanifold edges do not alone prove printability. MMGS uses default feature detection and no prescribed matching of opposite caps.</p>
<p><a href="surface_tradeoffs.json">Raw measurements and calibration</a> · <a href="surface_tradeoffs.py">Benchmark</a> · <a href="surface_tradeoffs_relaxed.py">Relaxed MMGS comparison</a> · <a href="surface_tradeoffs_api.py">Microgen API timing</a></p></html>"""
)
page = page.replace(
    "<h2>Choose the work budget</h2>",
    """<h2>What the measurements show</h2>
<p>Meshers quality mode is faster than both MMGS settings on all four tested cases. Relaxed MMGS gives gyroids better first-percentile angles and more uniform sizes; meshers improves Split-P angles and gives lower sampled distance error than relaxed MMGS on all four cases. Tight MMGS gives slightly lower sampled distance error on gyroids, at much higher counts and cost. These are different tradeoffs, not a universal quality ranking.</p>
<p>All outputs are closed with zero nonmanifold edges. Meshers matches the requested periodic faces. MMGS changes gyroid cap nodes and triangles. MMGS also changes Split-P cap nodes and triangles when checked on the original requested box planes. No unused vertices explain those mismatches. Neither MMGS setting passes the complete periodic check.</p>
<p>Fast meshers beats raw VTK in all four cases, with worse sampled geometric error. Relaxed MMGS counts are close for gyroids, but remain 25–32% above VTK for Split-P; meshers quality counts range from 10% below to 7% above VTK. The comparison is not exactly count-matched. Timings show medians with the three-run minimum and maximum beneath them; later API and relaxed-MMGS blocks were not globally interleaved, so small differences are not decisive.</p>
<h2>Choose the work budget</h2>""",
)
(root / "surface_tradeoffs.html").write_text(page, encoding="utf-8")
print(root / "surface_tradeoffs.html")
