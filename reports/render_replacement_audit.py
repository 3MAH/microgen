"""Build the replacement audit from measured data and the gap registry."""

import html
import json
import math
from pathlib import Path

root = Path(__file__).resolve().parent


def load(name):
    return json.loads((root / name).read_text())


gaps = load("replacement_gaps.json")
fidelity = load("replacement_fidelity.json")
fea = load("replacement_fedoo.json")
success = sorted(
    [r for r in fea if r["status"] == "ok" and r["mode"] == "meshers"],
    key=lambda r: r["resolution"],
)
intro = """The current experiment does not meet the replacement goal. Speed gains coexist with geometric-fidelity losses, inferior angle/size distributions in some cases, incomplete integration and unverified FEA accuracy. This audit inventories the known gaps across the current integration and marks untested paths explicitly. It does not claim every possible geometry and parameter combination has been benchmarked."""
contract = """Primary comparisons must use the same physical geometry and final triangle budget, or the same tetrahedron/DOF budget for solids. Counts and reference uncertainty must be reported, and any tolerance must be explicit. Surface quality requires angle and size distributions AND geometric fidelity. Simulation quality requires actual fedoo convergence, reactions, energy and solver checks. Fast/printing surfaces require fidelity, wall thickness and topology, even when element shape is unimportant. A failed baseline is a failure record, not a speed or accuracy win."""
sequence = [
    (
        "1. Fix the comparison contract",
        "Validate baseline topology and patch tests; implement explicit count budgets; retain failures and all costs. Add error-versus-count and time-to-error curves. Repair the MMG affine-patch baseline before using it for FEA comparisons.",
    ),
    (
        "2. Improve fast geometry",
        "Test bounded one/two-step edge-root corrections and error-triggered refinement. Measure facet interiors and reverse distances. If the six-tetrahedron decomposition wastes the budget, prototype direct hexahedral extraction with consistent ambiguity handling and periodic caps.",
    ),
    (
        "3. Improve surface shape at fixed count",
        "Implement periodic-orbit paired flips and tangential relocation, then count-neutral split/collapse redistribution. Optimize size distribution and angle tails together while rejecting geometric-error regressions. Start with the unit gyroid and graded Split-P tail defect.",
    ),
    (
        "4. Validate actual mechanics",
        "For solids: affine patch, manufactured solution, compression, then six periodic homogenization loads on a refinement series. For shells: first generate genuine midsurfaces, then bending/membrane cases with the same thickness, formulation and BCs. Use a converged reference and solve-cost accounting.",
    ),
    (
        "5. Improve scaling and fill integration gaps",
        "Use work-based worker selection, native field evaluation, deterministic block extraction, rolling caches and memory budgets. Then close anisotropy, chart seams/poles, sweeps, infill, density fitting, part types and generic-shape coverage. Release only the paths that meet their acceptance criteria.",
    ),
]
md = [
    "# Meshers replacement audit",
    intro,
    "## Acceptance contract",
    contract,
    "## Implementation order",
]
md += [f"### {a}\n\n{b}" for a, b in sequence]
md += ["## Complete gap registry"]
for g in gaps:
    md.append(
        f"### {g['id']} · {g['area']} · {g['priority']} · {g['status']}\n\n{g['problem']}\n\nEvidence: {g['evidence']}\n\nImprove: {g['improvement']}\n\nAcceptance: {g['acceptance']}"
    )
md += [
    "## Fedoo findings",
    "Five meshers gyroid refinement levels pass affine patch and manufactured-solution checks. Compression stiffness is not yet converged: the final 28-to-36 sample refinement changes stiffness by about 3.3%. The 16-sample stiffness is about 15% above the finest tested value, which is itself not an exact reference. This is clamped compression with free side surfaces, not periodic homogenization. No shell-FEA result is claimed.",
    "The 16-sample MMG baseline failed both fedoo and independent affine tests; four internal faces have same-side adjacent tetrahedra. Double-precision reruns did not fix it. The 20-sample fedoo patch also failed; its topology has not yet been diagnosed. The 12-sample attempt timed out after 180 seconds. These outputs do not provide a valid comparative FEA reference.",
]
(root / "meshers_replacement_audit.md").write_text("\n\n".join(md), encoding="utf-8")
escape = html.escape
cards = "".join(
    f'''<article data-status="{escape(g["status"])}" data-priority="{g["priority"]}"><div class="tag">{g["priority"]} · {escape(g["status"])}</div><h3>{g["id"]}. {escape(g["area"])}</h3><p>{escape(g["problem"])}</p><details><summary>Evidence, change and acceptance test</summary><p><b>Evidence.</b> {escape(g["evidence"])}</p><p><b>Change.</b> {escape(g["improvement"])}</p><p><b>Acceptance.</b> {escape(g["acceptance"])}</p></details></article>'''
    for g in gaps
)
names = {
    "vtk": "VTK",
    "fast": "Meshers fast",
    "accurate": "Meshers accurate",
    "quality": "Meshers quality",
    "mmgs_count": "VTK + MMGS relaxed",
}
colors = {
    "vtk": "#9aabba",
    "fast": "#57d6ca",
    "accurate": "#74a1ff",
    "quality": "#c59aff",
    "mmgs_count": "#ffbd70",
}
plots = []
tables = []
for case in dict.fromkeys(r["case"] for r in fidelity):
    rows = [r for r in fidelity if r["case"] == case]
    scale = max(r["p95"] for r in rows)
    svg = [
        '<svg viewBox="0 0 700 225" role="img" aria-label="Sampled surface distance by workflow">'
    ]
    for i, r in enumerate(rows):
        width = r["p95"] / scale * 390
        y = 12 + i * 39
        svg.append(
            f'<text x="0" y="{y + 17}">{names[r["mode"]]}</text><rect x="200" y="{y}" width="{width}" height="25" rx="4" fill="{colors[r["mode"]]}"/><text x="{208 + width}" y="{y + 17}">{r["p95"]:.5f}</text>'
        )
    svg.append("</svg>")
    body = "".join(
        f"<tr><th>{names[r['mode']]}</th><td>{r['triangles']:,}</td><td>{r['p95']:.5f}</td><td>{r['p99']:.5f}</td><td>{r['max_sampled']:.5f}</td></tr>"
        for r in rows
    )
    check = rows[0]["reference_check"]
    plots.append(
        f'<section><h3>{case.replace("_", " ")}</h3>{"".join(svg)}<p>Bidirectional p95 surface distance. Lower is better.</p><table><tr><th>Path</th><th>Triangles</th><th>p95</th><th>p99</th><th>Worst sampled</th></tr>{body}</table><p class="muted">Reference refinement p95 {check["p95"]:.5f}, worst sampled {check["max_sampled"]:.5f}. Counts {check["triangles"][0]:,} and {check["triangles"][1]:,}.</p></section>'
    )
fea_rows = "".join(
    f"<tr><td>{r['resolution']}</td><td>{r['tetrahedra']:,}</td><td>{r['minimum_quality']:.3f}</td><td>{r['mms']['l2_displacement_error']:.6f}</td><td>{r['mms']['energy_error']:.5f}</td><td>{r['compression']['stiffness']:.6f}</td><td>{r['compression']['cg_iterations']}</td></tr>"
    for r in success
)
xy = [
    (
        70
        + math.log(r["tetrahedra"] / success[0]["tetrahedra"])
        / math.log(success[-1]["tetrahedra"] / success[0]["tetrahedra"])
        * 550,
        210 - (r["compression"]["stiffness"] - 0.12) / 0.045 * 160,
    )
    for r in success
]
line = " ".join(f"{x:.1f},{y:.1f}" for x, y in xy)
fea_plot = f'<svg viewBox="0 0 700 270" role="img" aria-label="Compression stiffness versus refinement"><polyline points="{line}" fill="none" stroke="#c59aff" stroke-width="3"/>'
for r, (x, y) in zip(success, xy):
    fea_plot += f'<circle cx="{x}" cy="{y}" r="5" fill="#c59aff"/><text x="{x - 25}" y="{y - 15}">{r["compression"]["stiffness"]:.4f}</text><text x="{x - 20}" y="245">{r["tetrahedra"]:,}</text>'
fea_plot += '<text x="260" y="265">Tetrahedra, logarithmic spacing</text></svg>'
order = "".join(f"<h3>{a}</h3><p>{b}</p>" for a, b in sequence)
page = f"""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Meshers replacement audit</title><style>body{{background:#101821;color:#e6edf5;font:16px/1.55 system-ui;max-width:1450px;margin:auto;padding:32px}}h1{{font-size:40px;line-height:1.15}}h2{{margin-top:40px}}a{{color:#7ddadd}}section,article{{background:#192634;border:1px solid #344b60;border-radius:12px;padding:22px;margin:20px 0}}article h3{{margin:8px 0}}.tag,.muted{{color:#a9bdd0;font-size:13px}}.note{{border-left:5px solid #f3b66c;padding:16px 24px;background:#293445}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(450px,1fr));gap:20px}}.grid section{{margin:0}}table{{width:100%;border-collapse:collapse;font-size:13px}}td,th{{text-align:left;padding:9px;border-bottom:1px solid #344b60}}svg{{width:100%;max-width:750px}}svg text{{fill:#c5d4e3;font-size:12px}}input,select{{background:#1b2c3e;color:#e6edf5;border:1px solid #49667e;padding:12px;border-radius:8px;font:inherit}}input{{width:min(65%,650px)}}summary{{cursor:pointer;color:#8dd7e8}}p{{max-width:1100px}}[hidden]{{display:none!important}}@media(max-width:650px){{body{{padding:16px}}.grid{{display:block}}table{{display:block;overflow:auto}}}}</style>
<h1>Meshers is not yet a complete improvement over microgen</h1><p>29 September 2026 · 32 known gaps and unverified requirements · experimental source build</p><p class="note">{intro}</p><h2>The replacement criterion</h2><p>{contract}</p><p>Exact-count superiority has not been demonstrated. The plots below retain actual counts and expose losses; they are diagnostic evidence, not acceptance certificates.</p>
<h2>Geometric fidelity: stronger checks reveal losses</h2><p>Each mesh is compared in both directions with a fine exact-edge reference using closest-triangle distances from 10,000 area-weighted samples per direction. Reference grids use 63 and 95 cells per domain axis, with all optimization disabled. Distances are in model units; one TPMS period is one unit. This samples full closed surfaces, including caps. It is not a certified Hausdorff bound or a wall-thickness test. References use the same mesher, so an independently generated reference is still required for release validation.</p><div class="grid">{"".join(plots)}</div>
<h2>Fedoo: convergence, not just triangle or tetrahedron scores</h2><p>Fedoo 1.0.0, tet4 elasticity, E=2.5 and nu=0.25. Manufactured-solution and affine patch checks prescribe analytic displacement on the complete boundary. The separate structural test clamps the bottom and prescribes 0.001 axial displacement at the top with free sides. It audits energy/reaction agreement and the free-equation residual. Five meshers levels pass the checks.</p>{fea_plot}<table><tr><th>Samples/axis</th><th>Tetrahedra</th><th>Minimum shape quality</th><th>MMS L2 error</th><th>MMS energy error</th><th>Compression stiffness</th><th>CG iterations</th></tr>{fea_rows}</table>
<p class="note">The final stiffness still changes by about 3.3% between the last two levels. The common 16-sample mesh is about 15% stiffer than the finest tested mesh. That finest value is not converged ground truth. Good element shape does not establish adequate structural accuracy.</p><p>The MMG 16-sample output fails the affine patch test in both fedoo and an independent assembly. Four interior faces have adjacent tetrahedra on the same side, indicating local overlap. Double precision does not resolve it. The 20-sample fedoo patch also fails; that case's cause remains uninvestigated. The 12-sample attempt timed out after 180 seconds. These failures prevent a valid comparative FEA verdict. No shell or periodic-homogenization validation of the new surfaces is claimed.</p>
<h2>Complete gap registry</h2><p>Measured losses, missing adapters and untested behavior are different statuses. A returned native mesh is not automatically a supported microgen workflow. Filter by priority or search geometry, failure, evidence or proposed change.</p><input id="search" placeholder="Search all 32 gaps" aria-label="Search gaps"><select id="priority" aria-label="Priority"><option value="">All priorities</option><option>P0</option><option>P1</option><option>P2</option></select><div id="gaps">{cards}</div>
<h2>Implementation order</h2>{order}
<h2>Evidence and reproduction</h2><p><a href="meshers_replacement_audit.md">Full text audit</a> · <a href="replacement_gaps.json">Gap registry</a> · <a href="replacement_fidelity.json">Distance data</a> · <a href="replacement_fedoo.json">Fedoo data</a> · <a href="replacement_baseline_audit.json">Baseline defect</a> · <a href="surface_tradeoffs.html">Earlier timings and triangle statistics</a></p><p>Scripts: replacement_fidelity.py, replacement_fedoo.py, replacement_fedoo_refine.py and check_fedoo_mmg.py. They use the experimental meshers package, its examples and demos/solver on PYTHONPATH, plus the MMGS/MMG3D binaries bundled with mmgpy 0.16.2 on PATH. Fedoo interfaces follow its <a href="https://3mah.github.io/fedoo-docs/Quick_Start.html">official documentation</a>. The new FEA timings are diagnostic runs under development load, not controlled speed comparisons.</p>
<script>const q=document.querySelector('#search'),p=document.querySelector('#priority');function filter(){{for(const a of document.querySelectorAll('#gaps article'))a.hidden=!(a.textContent.toLowerCase().includes(q.value.toLowerCase())&&(!p.value||a.dataset.priority===p.value));}}q.addEventListener('input',filter);p.addEventListener('change',filter);</script></html>"""
(root / "meshers_replacement_audit.html").write_text(page, encoding="utf-8")
print(root / "meshers_replacement_audit.html")
