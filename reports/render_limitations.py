"""Render the self-contained HTML report from limitations_data.json."""

import html
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "limitations_data.json").read_text(encoding="utf-8"))
CORE = json.loads((HERE / "curved_core_probe_data.json").read_text(encoding="utf-8"))
MATCHED = json.loads((HERE / "matched_surface_data.json").read_text(encoding="utf-8"))
FAST = json.loads((HERE / "fast_surface_data.json").read_text(encoding="utf-8"))
LINEAR = json.loads((HERE / "linear_surface_data.json").read_text(encoding="utf-8"))
TOPOLOGY = json.loads((HERE / "fast_surface_topology.json").read_text(encoding="utf-8"))
FAST_VOLUME = json.loads((HERE / "fast_volume_data.json").read_text(encoding="utf-8"))
RUNS = DATA["runs"]
SURFACE_CASES = [
    ("gyroid_unit", "Gyroid · 1 cell"),
    ("split_p_unit", "Split-P · 1 cell"),
    ("graded_gyroid_2", "Graded gyroid · 2³"),
    ("graded_gyroid_3", "Graded gyroid · 3³"),
    ("graded_split_p_2", "Graded split-P · 2³"),
    ("lateral_gyroid_2", "X-graded gyroid · 2³"),
]
VOLUME_CASES = [
    ("gyroid_unit", "Gyroid · 1 cell"),
    ("split_p_unit", "Split-P · 1 cell"),
    ("graded_gyroid_2", "Graded gyroid · 2³"),
    ("lateral_gyroid_2", "X-graded gyroid · 2³"),
]
CURVED_CASES = [
    ("cylinder_full_wrap", "Full-wrap cylinder"),
    ("sphere_full_wrap", "Full sphere"),
    ("sweep", "Sweep"),
]


def samples(group: str, case: str, mode: str) -> list[dict]:
    return [
        run
        for run in RUNS
        if run["group"] == group
        and run["case"] == case
        and run["mode"] == mode
        and run["status"] == "ok"
    ]


def median(group: str, case: str, mode: str, key: str = "runtime_seconds") -> float:
    return statistics.median(run[key] for run in samples(group, case, mode))


def representative(group: str, case: str, mode: str) -> dict:
    return samples(group, case, mode)[0]


def fmt(value: float) -> str:
    return f"{value:.3f}" if value < 10 else f"{value:.2f}"


def residual_fmt(value: float) -> str:
    return f"{value:.1e}" if value < 0.001 else f"{value:.3f}"


def matched_samples(case: str, mode: str) -> list[dict]:
    return [
        run
        for run in MATCHED["runs"]
        if run["phase"] == "timed"
        and run["case"] == case
        and run["mode"] == mode
        and run["status"] == "ok"
    ]


def matched_result(case: str, mode: str) -> dict:
    runs = matched_samples(case, mode)
    assert len(runs) == 3, (case, mode)
    return runs[0]


def matched_time(case: str, mode: str) -> float:
    return statistics.median(
        run["generation_seconds"] for run in matched_samples(case, mode)
    )


def speed_runs(data: dict, case: str, mode: str) -> list[dict]:
    runs = [
        run
        for run in data["runs"]
        if run["phase"] == "timed"
        and run["case"] == case
        and run["mode"] == mode
        and run["status"] == "ok"
    ]
    assert len(runs) == 3, (case, mode)
    return runs


def speed_time(data: dict, case: str, mode: str) -> float:
    return statistics.median(
        run["generation_seconds"] for run in speed_runs(data, case, mode)
    )


def fast_surface_table() -> str:
    lines = [
        '<div class="table-wrap"><table><thead><tr><th>Case</th>'
        "<th>Triangles VTK / fast meshers</th><th>VTK / fast exact / linear time</th>"
        "<th>Fast exact / VTK</th><th>Linear / VTK</th>"
        "<th>1st percentile min angle VTK / fast / linear</th>"
        "<th>Area CV VTK / fast / linear</th>"
        "<th>95th percentile field residual VTK / fast / linear</th></tr></thead><tbody>"
    ]
    for case, name in SURFACE_CASES:
        vtk = speed_runs(LINEAR, case, "vtk")[0]
        fast = speed_runs(FAST, case, "meshers_fast")[0]
        linear = speed_runs(LINEAR, case, "meshers_linear")[0]
        vtk_time = statistics.median(
            run["generation_seconds"]
            for data in (FAST, LINEAR)
            for run in speed_runs(data, case, "vtk")
        )
        fast_time = speed_time(FAST, case, "meshers_fast")
        linear_time = speed_time(LINEAR, case, "meshers_linear")
        lines.append(
            f"<tr><th>{html.escape(name)}</th>"
            f"<td>{vtk['triangle_count']:,} / {fast['triangle_count']:,}</td>"
            f"<td>{fmt(vtk_time)} / {fmt(fast_time)} / {fmt(linear_time)} s</td>"
            f"<td>{fast_time / vtk_time:.2f}×</td>"
            f"<td>{linear_time / vtk_time:.2f}×</td>"
            f"<td>{vtk['angle_p01_degrees']:.1f}° / {fast['angle_p01_degrees']:.1f}° / {linear['angle_p01_degrees']:.1f}°</td>"
            f"<td>{vtk['area_cv']:.2f} / {fast['area_cv']:.2f} / {linear['area_cv']:.2f}</td>"
            f"<td>{residual_fmt(vtk['implicit_vertex_residual_p95'])} / {residual_fmt(fast['implicit_vertex_residual_p95'])} / {residual_fmt(linear['implicit_vertex_residual_p95'])}</td></tr>"
        )
    lines.append("</tbody></table></div>")
    return "".join(lines)


def fast_speed_chart() -> str:
    width, height, left, right = 1060, 95 + 59 * len(SURFACE_CASES), 245, 945

    def xpos(ratio: float) -> float:
        return left + ratio / 2.5 * (right - left)

    parts = [
        f'<svg class="plot" viewBox="0 0 {width} {height}" role="img" '
        'aria-label="Unpolished and linear meshers generation time relative to raw VTK at similar triangle counts">',
        '<text x="245" y="23" class="axis-title">Generation time / raw VTK time · lower is faster</text>',
    ]
    for tick in (0, 0.5, 1, 1.5, 2, 2.5):
        x = xpos(tick)
        dash = ' stroke-dasharray="5 5" stroke="#eaf3f3"' if tick == 1 else ""
        parts.append(
            f'<line x1="{x:.1f}" x2="{x:.1f}" y1="40" y2="{height - 34}" class="grid"{dash}/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{height - 10}" text-anchor="middle" class="tick">{tick:g}×</text>'
        )
    for index, (case, name) in enumerate(SURFACE_CASES):
        y = 68 + index * 59
        vtk = statistics.median(
            run["generation_seconds"]
            for data in (FAST, LINEAR)
            for run in speed_runs(data, case, "vtk")
        )
        fast = speed_time(FAST, case, "meshers_fast") / vtk
        linear = speed_time(LINEAR, case, "meshers_linear") / vtk
        parts.append(
            f'<text x="10" y="{y + 4}" class="row-label">{html.escape(name)}</text>'
        )
        parts.append(
            f'<line x1="{xpos(fast):.1f}" x2="{xpos(linear):.1f}" y1="{y}" y2="{y}" class="pair-line"/>'
        )
        parts.append(
            f'<circle cx="{xpos(fast):.1f}" cy="{y}" r="7" class="meshers-dot"><title>Exact intersections: {fast:.2f}× VTK</title></circle>'
        )
        parts.append(
            f'<circle cx="{xpos(linear):.1f}" cy="{y}" r="7" fill="var(--blue)"><title>Linear intersections: {linear:.2f}× VTK</title></circle>'
        )
    parts.append("</svg>")
    return "".join(parts)


def fast_volume_table() -> str:
    lines = [
        '<div class="table-wrap"><table><thead><tr><th>Case</th>'
        "<th>Raw VTK mixed cells</th><th>Meshers tetrahedra</th>"
        "<th>Time VTK / meshers</th><th>Meshers / VTK</th>"
        "<th>Meshers min MMG quality, optimized / unoptimized</th></tr></thead><tbody>"
    ]
    for case, name in VOLUME_CASES:
        if case not in ("gyroid_unit", "split_p_unit", "graded_gyroid_2"):
            continue
        vtk = [
            run
            for run in FAST_VOLUME["runs"]
            if run["case"] == case and run["mode"] == "microgen_volume"
        ]
        meshers = [
            run
            for run in FAST_VOLUME["runs"]
            if run["case"] == case and run["mode"] == "meshers_volume"
        ]
        assert len(vtk) == len(meshers) == 3
        a = statistics.median(run["runtime_seconds"] for run in vtk)
        b = statistics.median(run["runtime_seconds"] for run in meshers)
        optimized = representative("volume", case, "meshers_volume")[
            "minimum_mmg_quality"
        ]
        unoptimized = meshers[0]["minimum_mmg_quality"]
        lines.append(
            f"<tr><th>{html.escape(name)}</th><td>{vtk[0]['elements']:,}</td>"
            f"<td>{meshers[0]['elements']:,}</td><td>{fmt(a)} / {fmt(b)} s</td>"
            f"<td>{b / a:.1f}×</td><td>{optimized:.3f} / {unoptimized:.3f}</td></tr>"
        )
    lines.append("</tbody></table></div>")
    return "".join(lines)


def matched_table() -> str:
    lines = [
        '<div class="table-wrap"><table><thead><tr><th>Case</th><th>Grid points/cell VTK / meshers</th><th>Triangles VTK / meshers</th>'
        "<th>Count gap</th><th>Time VTK / meshers</th><th>Time ratio</th>"
        "<th>1st percentile min angle VTK / meshers</th>"
        "<th>Triangles below 10° VTK / meshers</th>"
        "<th>Area CV VTK / meshers</th></tr></thead><tbody>"
    ]
    for case, name in SURFACE_CASES:
        a, b = matched_result(case, "vtk"), matched_result(case, "meshers")
        ta, tb = matched_time(case, "vtk"), matched_time(case, "meshers")
        gap = 100 * (b["triangle_count"] / a["triangle_count"] - 1)
        lines.append(
            f"<tr><th>{html.escape(name)}</th>"
            f"<td>{a['resolution']} / {b['resolution']}</td>"
            f"<td>{a['triangle_count']:,} / {b['triangle_count']:,}</td>"
            f"<td>{gap:+.1f}%</td><td>{fmt(ta)} / {fmt(tb)} s</td>"
            f"<td>{tb / ta:.1f}×</td>"
            f"<td>{a['angle_p01_degrees']:.1f}° / {b['angle_p01_degrees']:.1f}°</td>"
            f"<td>{100 * a['fraction_min_angle_below_10_degrees']:.1f}% / {100 * b['fraction_min_angle_below_10_degrees']:.1f}%</td>"
            f"<td>{a['area_cv']:.2f} / {b['area_cv']:.2f}</td></tr>"
        )
    lines.append("</tbody></table></div>")
    return "".join(lines)


def matched_size_chart() -> str:
    # Equivalent-area diameter, divided by the VTK median in each case.
    width, height, left, right = 1060, 85 + 82 * len(SURFACE_CASES), 245, 950
    maximum = max(
        matched_result(case, mode)["equivalent_diameter_p99"]
        / matched_result(case, "vtk")["equivalent_diameter_median"]
        for case, _ in SURFACE_CASES
        for mode in ("vtk", "meshers")
    )
    xhi = max(2.0, math.ceil(maximum * 2) / 2)

    def xpos(value):
        return left + value / xhi * (right - left)

    parts = [
        f'<svg class="plot" viewBox="0 0 {width} {height}" role="img" '
        'aria-label="Triangle equivalent diameter distributions; first percentile, median and ninety-ninth percentile relative to each VTK median">',
        '<text x="245" y="22" class="axis-title">Equivalent-area diameter / VTK median within each case</text>',
    ]
    for tick in (0, 0.5, 1, 1.5, 2, 2.5, 3):
        if tick > xhi:
            break
        x = xpos(tick)
        parts.append(
            f'<line x1="{x:.1f}" x2="{x:.1f}" y1="38" y2="{height - 32}" class="grid"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{height - 9}" text-anchor="middle" class="tick">{tick:g}</text>'
        )
    for index, (case, name) in enumerate(SURFACE_CASES):
        y = 69 + index * 82
        reference = matched_result(case, "vtk")["equivalent_diameter_median"]
        parts.append(
            f'<text x="10" y="{y + 13}" class="row-label">{html.escape(name)}</text>'
        )
        for offset, mode, color in (
            (0, "vtk", "var(--vtk)"),
            (24, "meshers", "var(--mesh)"),
        ):
            run = matched_result(case, mode)
            lo, mid, hi = (
                xpos(run[f"equivalent_diameter_{key}"] / reference)
                for key in ("p01", "median", "p99")
            )
            yy = y + offset
            parts.append(
                f'<line x1="{lo:.1f}" x2="{hi:.1f}" y1="{yy}" y2="{yy}" stroke="{color}" stroke-width="5"/>'
            )
            parts.append(f'<circle cx="{mid:.1f}" cy="{yy}" r="6" fill="{color}"/>')
    parts.append("</svg>")
    return "".join(parts)


def log_chart(group: str, cases: list[tuple[str, str]], label: str) -> str:
    width, height = 1060, 100 + 60 * len(cases)
    left, right = 250, 930
    xlo, xhi = -2.0, 1.0

    def x(value: float) -> float:
        return left + (math.log10(value) - xlo) / (xhi - xlo) * (right - left)

    parts = [
        f'<svg class="plot" viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="{html.escape(label)}; time in seconds on a logarithmic scale">',
        f'<text x="{left}" y="22" class="axis-title">Generation time · seconds · log scale</text>',
    ]
    for tick in (0.01, 0.1, 1, 10):
        xx = x(tick)
        parts.append(
            f'<line x1="{xx:.1f}" x2="{xx:.1f}" y1="42" y2="{height - 32}" class="grid"/>'
        )
        parts.append(
            f'<text x="{xx:.1f}" y="{height - 12}" text-anchor="middle" class="tick">{tick:g}</text>'
        )
    for i, (case, name) in enumerate(cases):
        yy = 68 + 60 * i
        a = median(group, case, f"microgen_{group}")
        b = median(group, case, f"meshers_{group}")
        ax, bx = x(a), x(b)
        parts.append(
            f'<text x="12" y="{yy + 4}" class="row-label">{html.escape(name)}</text>'
        )
        parts.append(
            f'<line x1="{ax:.1f}" x2="{bx:.1f}" y1="{yy}" y2="{yy}" class="pair-line"/>'
        )
        parts.append(
            f'<circle cx="{ax:.1f}" cy="{yy}" r="8" class="vtk-dot"><title>microgen VTK: {fmt(a)} s</title></circle>'
        )
        parts.append(
            f'<circle cx="{bx:.1f}" cy="{yy}" r="8" class="meshers-dot"><title>meshers: {fmt(b)} s</title></circle>'
        )
        parts.append(
            f'<text x="{right + 18}" y="{yy + 5}" class="ratio">{b / a:.1f}×</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


def comparison_table(group: str, cases: list[tuple[str, str]]) -> str:
    lines = [
        '<div class="table-wrap"><table><thead><tr><th>Case</th><th>microgen VTK</th>'
        "<th>meshers</th><th>meshers / VTK</th><th>Output</th></tr></thead><tbody>"
    ]
    for case, name in cases:
        a = median(group, case, f"microgen_{group}")
        b = median(group, case, f"meshers_{group}")
        va = representative(group, case, f"microgen_{group}")
        vb = representative(group, case, f"meshers_{group}")
        if group == "surface":
            output = f"{va['elements']:,} vs {vb['elements']:,} triangles"
        else:
            output = f"{va['elements']:,} mixed cells vs {vb['elements']:,} tets"
        lines.append(
            f"<tr><th>{html.escape(name)}</th><td>{fmt(a)} s</td><td>{fmt(b)} s</td>"
            f"<td>{b / a:.1f}×</td><td>{output}</td></tr>"
        )
    lines.append("</tbody></table></div>")
    return "".join(lines)


def curved_table() -> str:
    lines = [
        '<div class="table-wrap"><table><thead><tr><th>Geometry</th><th>VTK surface</th>'
        "<th>VTK volume</th><th>Current adapter</th><th>Surface triangles</th>"
        "<th>Open edges</th></tr></thead><tbody>"
    ]
    for case, name in CURVED_CASES:
        surface = median("curved", case, "microgen_surface")
        volume = median("curved", case, "microgen_volume")
        surface_result = representative("curved", case, "microgen_surface")
        count = surface_result["elements"]
        lines.append(
            f"<tr><th>{html.escape(name)}</th><td>{fmt(surface)} s</td>"
            f'<td>{fmt(volume)} s</td><td><span class="tag">VTK fallback</span></td>'
            f"<td>{count:,}</td><td>{surface_result['open_edges']:,}</td></tr>"
        )
    lines.append("</tbody></table></div>")
    return "".join(lines)


def core_probe(case: str, kind: str) -> dict:
    return next(
        run for run in CORE["runs"] if run["case"] == case and run["kind"] == kind
    )


def core_table() -> str:
    lines = [
        '<div class="table-wrap"><table><thead><tr><th>Geometry</th><th>Meshers volume probe</th><th>Meshers surface probe</th></tr></thead><tbody>'
    ]
    for case, name in (
        ("cylinder_sector", "Cylinder sector"),
        ("sphere_sector", "Sphere sector"),
        *CURVED_CASES,
    ):
        volume = core_probe(case, "volume")
        surface = core_probe(case, "surface")
        if volume["status"] == "mesh_returned":
            volume_text = (
                f"{volume['tetrahedra']:,} tets · {fmt(volume['seconds'])} s · "
                f"qmin {volume['minimum_mmg_quality']:.3f} · "
                f"error {volume['sampled_surface_error']:.3f}"
            )
            if case == "cylinder_full_wrap":
                volume_text += (
                    f" · welded open edges {volume['welded_volume_surface_open_edges']}"
                )
        else:
            volume_text = f"Rejected: {volume['reason']}"
        if surface["status"] == "mesh_returned":
            surface_text = (
                f"{surface['triangles']:,} triangles · {fmt(surface['seconds'])} s · "
                f"post-map open edges {surface['open_edges_after_seam_cleanup']:,}"
            )
        else:
            surface_text = f"Rejected: {surface['reason']}"
        lines.append(
            f"<tr><th>{html.escape(name)}</th><td>{html.escape(volume_text)}</td>"
            f"<td>{html.escape(surface_text)}</td></tr>"
        )
    lines.append("</tbody></table></div>")
    return "".join(lines)


def raw_times(group: str, cases: list[tuple[str, str]]) -> str:
    lines = []
    for case, name in cases:
        modes = (
            ("microgen_surface", "microgen_volume")
            if group == "curved"
            else (f"microgen_{group}", f"meshers_{group}")
        )
        for mode in modes:
            values = ", ".join(
                fmt(run["runtime_seconds"]) for run in samples(group, case, mode)
            )
            lines.append(
                f"<tr><th>{html.escape(name)}</th><td>{mode}</td><td>{values}</td></tr>"
            )
    return "".join(lines)


def validate_data() -> None:
    assert (
        FAST["meshers_commit"]
        == LINEAR["meshers_commit"]
        == FAST_VOLUME["meshers_commit"]
        == TOPOLOGY["meshers_commit"]
    )
    assert len(TOPOLOGY["runs"]) == 2 * len(SURFACE_CASES)
    assert all(
        run["open_edges"] == 0
        and not any(run["periodic_node_mismatch"])
        and not any(run["periodic_triangle_mismatch"])
        for run in TOPOLOGY["runs"]
    )
    for group, cases in (("surface", SURFACE_CASES), ("volume", VOLUME_CASES)):
        for case, _ in cases:
            for mode in (f"microgen_{group}", f"meshers_{group}"):
                assert len(samples(group, case, mode)) == 3, (group, case, mode)
    for case, _ in CURVED_CASES:
        for mode in ("microgen_surface", "microgen_volume"):
            assert len(samples("curved", case, mode)) == 3, (case, mode)
        unsupported = [
            run
            for run in RUNS
            if run["group"] == "curved"
            and run["case"] == case
            and run["mode"] == "meshers_volume"
        ]
        assert len(unsupported) == 1 and unsupported[0]["status"] == "unsupported"


def build() -> str:
    validate_data()
    surface_chart = log_chart("surface", SURFACE_CASES, "Surface generation comparison")
    volume_chart = log_chart("volume", VOLUME_CASES, "Raw volume generation comparison")
    surface_table = comparison_table("surface", SURFACE_CASES)
    volume_table = comparison_table("volume", VOLUME_CASES)
    low_quality_volume_table = fast_volume_table()
    coverage_table = curved_table()
    capability_table = core_table()
    matched_comparison_table = matched_table()
    fast_comparison_table = fast_surface_table()
    fast_chart = fast_speed_chart()
    size_chart = matched_size_chart()
    matched_gaps = [
        abs(
            matched_result(case, "meshers")["triangle_count"]
            / matched_result(case, "vtk")["triangle_count"]
            - 1
        )
        for case, _ in SURFACE_CASES
    ]
    polished_ratios = [
        matched_time(case, "meshers") / matched_time(case, "vtk")
        for case, _ in SURFACE_CASES
    ]
    fast_wins = sum(
        speed_time(FAST, case, "meshers_fast")
        < statistics.median(
            run["generation_seconds"]
            for data in (FAST, LINEAR)
            for run in speed_runs(data, case, "vtk")
        )
        for case, _ in SURFACE_CASES
    )
    sample_rows = raw_times("surface", SURFACE_CASES) + raw_times(
        "volume", VOLUME_CASES
    )
    date = html.escape(DATA["measured_at_utc"][:10])
    microgen_rev = html.escape(DATA["microgen_commit"][:10])
    meshers_rev = html.escape(DATA["meshers_commit"][:10])
    prototype_rev = html.escape(FAST["meshers_commit"][:10])
    platform_name = html.escape(DATA["platform"])
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>microgen × meshers | Three limits to a full switch</title>
<style>
:root {{ --bg:#0c1420; --panel:#152333; --panel2:#1b2d40; --ink:#eaf3f3; --muted:#a7bac5; --line:#365064; --vtk:#65d6c4; --mesh:#ffae79; --blue:#80aaff; }}
* {{ box-sizing:border-box }} html {{ scroll-behavior:smooth }} body {{ margin:0; background:var(--bg); color:var(--ink); font:16px/1.58 system-ui,-apple-system,Segoe UI,sans-serif }}
a {{ color:#9cd7ff }} a:hover {{ color:white }} .shell {{ max-width:1160px; margin:auto; padding:0 30px 90px }}
header {{ padding:75px 0 40px; border-bottom:1px solid var(--line) }} .eyebrow {{ text-transform:uppercase; letter-spacing:.18em; color:var(--vtk); font-size:.78rem; font-weight:800 }}
h1 {{ font-size:clamp(2.4rem,5vw,4.5rem); line-height:1.07; letter-spacing:-.045em; max-width:880px; margin:16px 0 22px }} h2 {{ font-size:clamp(1.7rem,3vw,2.5rem); letter-spacing:-.03em; line-height:1.15; margin:0 0 20px }} h3 {{ font-size:1.15rem; margin:0 0 8px }}
.lead {{ font-size:1.23rem; color:#c7d9dc; max-width:850px; margin:0 0 28px }} .meta {{ display:flex; flex-wrap:wrap; gap:14px; color:var(--muted); font-size:.9rem }} .meta span {{ border:1px solid var(--line); border-radius:40px; padding:5px 13px }}
nav {{ display:flex; gap:12px; flex-wrap:wrap; padding:20px 0 }} nav a {{ text-decoration:none; padding:7px 12px; border-radius:8px; background:var(--panel) }}
.takeaway {{ margin:18px 0 62px; padding:25px 29px; background:linear-gradient(115deg,#173e45,#183149 62%,#2a3040); border:1px solid #3f6770; border-radius:18px }} .takeaway strong {{ color:#a7fff1 }}
.cards {{ display:grid; grid-template-columns:repeat(3,1fr); gap:15px; margin:22px 0 8px }} .card {{ background:var(--panel); border:1px solid var(--line); border-radius:14px; padding:20px }} .card .number {{ font-size:2.2rem; font-weight:800; letter-spacing:-.04em; color:var(--vtk); line-height:1.1 }} .card small {{ color:var(--muted); display:block; margin-top:8px }}
section {{ margin:75px 0 0 }} .section-tag {{ font-size:.8rem; text-transform:uppercase; color:var(--vtk); letter-spacing:.16em; font-weight:800; margin-bottom:12px }} p {{ max-width:900px }} .muted {{ color:var(--muted) }}
.chart-card {{ background:var(--panel); border:1px solid var(--line); border-radius:17px; margin:26px 0 18px; padding:18px 17px 12px }} .plot {{ width:100%; height:auto; display:block }} .grid {{ stroke:#334a59; stroke-width:1 }} .axis-title,.tick {{ fill:var(--muted); font:13px system-ui,sans-serif }} .row-label {{ fill:var(--ink); font:15px system-ui,sans-serif }} .ratio {{ fill:var(--ink); font:700 16px system-ui,sans-serif }} .pair-line {{ stroke:#6b8192; stroke-width:3 }} .vtk-dot {{ fill:var(--vtk); stroke:#0e2630; stroke-width:2 }} .meshers-dot {{ fill:var(--mesh); stroke:#30251f; stroke-width:2 }}
.legend {{ display:flex; flex-wrap:wrap; gap:22px; color:var(--muted); font-size:.9rem; margin:0 8px 8px }} .swatch {{ display:inline-block; width:11px; height:11px; border-radius:50%; margin-right:7px }}
.table-wrap {{ overflow-x:auto; border:1px solid var(--line); border-radius:12px }} table {{ border-collapse:collapse; width:100%; min-width:680px }} th,td {{ padding:11px 14px; text-align:left; border-bottom:1px solid #30465a; white-space:nowrap }} thead {{ color:var(--muted); font-size:.79rem; text-transform:uppercase; letter-spacing:.05em; background:#203144 }} tbody tr:last-child th,tbody tr:last-child td {{ border-bottom:0 }} tbody th {{ font-weight:600 }} tbody tr:nth-child(even) {{ background:#172838 }}
.tag {{ display:inline-block; background:#543a35; color:#ffd0aa; padding:3px 9px; border-radius:20px; font-size:.84rem; font-weight:700 }}
.note {{ border-left:3px solid var(--blue); background:#172738; padding:14px 18px; margin:20px 0; border-radius:0 9px 9px 0 }}
.flow {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; margin-top:24px }} .flow > div {{ border:1px solid var(--line); border-radius:12px; padding:18px; background:var(--panel) }} .steps {{ color:var(--muted); font-size:.92rem }} .steps b {{ color:var(--ink) }} .arrow {{ color:var(--vtk); padding:0 9px }}
details {{ background:var(--panel); border:1px solid var(--line); border-radius:12px; padding:16px 20px; margin:24px 0 }} summary {{ cursor:pointer; font-weight:700 }} details table {{ margin-top:15px }} code {{ font-family:ui-monospace,Consolas,monospace; font-size:.89em; color:#a8dcff }}
footer {{ margin-top:80px; padding-top:24px; border-top:1px solid var(--line); color:var(--muted) }}
@media(max-width:720px) {{ .shell {{ padding:0 16px 60px }} header {{ padding-top:42px }} .cards,.flow {{ grid-template-columns:1fr }} .chart-card {{ overflow-x:auto }} .plot {{ min-width:720px }} }}
</style>
</head>
<body><div class="shell">
<header>
<div class="eyebrow">Experimental benchmark report · {date}</div>
<h1>Three limits to replacing microgen’s mesh paths with meshers</h1>
<p class="lead">Quality-checked volume meshing should compare meshers with microgen’s VTK plus MMG workflow. This surface study asks a different question: can meshers beat raw VTK by relaxing triangle improvement? At similar triangle counts, disabling improvement wins only on the smallest gyroid. A linear-intersection prototype can get closer to VTK speed, but its geometric error rises sharply.</p>
<div class="meta"><span>Windows · Python {html.escape(DATA["python"])}</span><span>VTK 16 grid points per cell; meshers tuned</span><span>3 fresh-process trials per timed path</span><span>Imports excluded</span></div>
</header>
<nav><a href="#fast">01 Speed vs quality</a><a href="#matched">02 Similar triangle counts</a><a href="#surface">03 Fixed sampling</a><a href="#coverage">04 Curved coverage</a><a href="#volume">05 Raw volume speed</a><a href="#method">Methods</a></nav>
<div class="takeaway"><strong>Practical call.</strong> Keep the direct meshers path for quality-checked tetrahedral volumes, where VTK plus MMG is the relevant alternative. For print surfaces, dropping meshers’ triangle improvement makes it much faster but does not generally beat VTK. Linear edge intersections are faster still, yet their geometric error is too large in these probes, especially for split-P. The direct meshers surface generator is experimental and absent from the published 0.1.0 wheel.</div>
<div class="cards">
<div class="card"><div class="number">{fast_wins} / {len(SURFACE_CASES)}</div><div>Cases where unpolished meshers beats raw VTK</div><small>Similar triangle counts; exact edge intersections kept</small></div>
<div class="card"><div class="number">1 / 3</div><div>Current curved adapter fallbacks with a verified closed meshers volume path</div><small>Full-wrap cylinder passes; sphere and Sweep fail at singular axes</small></div>
<div class="card"><div class="number">{min(polished_ratios):.1f}–{max(polished_ratios):.1f}×</div><div>Polished meshers surface time over VTK at similar counts</div><small>Higher triangle quality in every tested case</small></div>
</div>

<section id="fast"><div class="section-tag">Surface experiment · speed versus quality</div><h2>Turning off improvement rarely beats raw VTK</h2>
<p>Meshers can skip all smoothing, local topology improvement and periodic polishing while keeping exact analytic edge intersections. That removes most of its surface-mesh work. The linear prototype also replaces exact edge roots with interpolation and reuses one normal per polygon. For each case, both fast modes were tuned independently to the nearest VTK triangle count. VTK used 16 points per cell. Each mode then ran three fresh-process trials; VTK time is the median of six trials across the two batches.</p>
{fast_comparison_table}
<div class="chart-card">{fast_chart}<div class="legend"><span><i class="swatch" style="background:var(--mesh)"></i>Unpolished exact intersections</span><span><i class="swatch" style="background:var(--blue)"></i>Linear intersections</span><span>Dashed threshold at 1×: equal to raw VTK time</span></div></div>
<div class="note">The unpolished exact-root mode beats VTK only for the one-cell gyroid. Linear intersections reduce time in some larger cases but still do not reliably beat VTK. They also lose geometric fidelity: the 95th-percentile residual of | |f| − 1 | at interior surface vertices is several times the VTK value and reaches 0.702 for unit split-P. This residual is dimensionless and compares the same normalized TPMS field. The fast modes mostly give up meshers’ angle and size-distribution advantage, although split-P retains some improvement. All 12 selected fast meshes had zero open edges, and their requested periodic caps matched. The linear path is an experiment, not a recommended microgen backend.</div>
<p class="muted">The exact-root fast mode retains near-zero residual at the measured interior vertices, while linear interpolation loses that accuracy. Vertex residual does not measure facet-interior deviation, wall thickness or printability.</p>
<p class="muted">One fast exact-root scan candidate failed at 10 points per cell for graded split-P with an inverted triangle; the selected 12-point setting passed all timed trials. Every selected linear setting also passed all timed trials.</p>
</section>

<section id="matched"><div class="section-tag">Surface comparison · matched output size</div><h2>Similar triangle counts reveal the real tradeoff</h2>
<p>For each case, VTK used 16 grid points per cell. The meshers input resolution was chosen from 10 to 15 points per cell to minimize the triangle-count gap. Each selected path then ran three times in a fresh process, with order alternated. Generation time excludes imports and the later metric calculations. The table reports medians for time and one deterministic mesh for shape and size metrics.</p>
{matched_comparison_table}
<p class="muted">Angle values are the 1st percentile of each triangle’s smallest interior angle. Area CV is standard deviation divided by mean, so lower means a tighter size distribution. Counts differ by at most {max(matched_gaps) * 100:.1f}%; this is a nearest-count comparison, not a promise of identical geometry error.</p>
<p class="muted">Two meshers scan candidates failed and were excluded from selection: graded split-P at 10 points per cell produced an inverted triangle; x-graded gyroid at 12 did not satisfy periodic cap matching and the 5° quality gate. All six selected settings succeeded in all three timed trials.</p>
<div class="chart-card">{size_chart}<div class="legend"><span><i class="swatch" style="background:var(--vtk)"></i>microgen VTK</span><span><i class="swatch" style="background:var(--mesh)"></i>meshers</span><span>Line = 1st to 99th percentile; dot = median</span></div></div>
<div class="note">The chart compares equivalent-area diameter, √(4A/π), after dividing by the VTK median within each case. A broad line means triangle areas span a broad range. These metrics describe triangle shape and size, not surface-position error, wall thickness or slicer behavior.</div>
</section>

<section id="surface"><div class="section-tag">Surface comparison · fixed sampling</div><h2>At equal input resolution, meshers emits more triangles</h2>
<p>microgen clips a 3D VTK grid and extracts its boundary. The meshers experiment extracts and improves triangles directly from a sampled 3D field. At matched grid-point counts, VTK generated fewer triangles and finished first in all six cases. The x-graded case retained matching y/z cap triangles in both outputs.</p>
<div class="chart-card">{surface_chart}<div class="legend"><span><i class="swatch" style="background:var(--vtk)"></i>microgen VTK</span><span><i class="swatch" style="background:var(--mesh)"></i>meshers direct surface</span><span>Right-hand number = meshers / VTK time</span></div></div>
{surface_table}
<div class="note">The VTK surfaces had zero open edges in these six tests. This checks closure, not printability: self-intersections, wall thickness, slicer behavior and dimensional error were not validated. Meshers’ higher minimum triangle angles are useful for numerical work, but this report does not treat angle as an additive-manufacturing acceptance criterion.</div>
<div class="flow"><div><h3>microgen VTK sheet</h3><div class="steps"><b>3D structured samples</b><span class="arrow">→</span><b>clip twice</b><span class="arrow">→</span><b>extract boundary triangles</b></div></div><div><h3>meshers experimental surface</h3><div class="steps"><b>3D sampled field</b><span class="arrow">→</span><b>surface triangles only</b><span class="arrow">→</span><b>improve / periodic polish</b></div></div></div>
</section>

<section id="coverage"><div class="section-tag">Limit 02 · geometry coverage</div><h2>The adapter fallback hides a working cylinder path</h2>
<p>microgen’s current adapter routes full cylindrical seams, spherical poles and Sweep to VTK. That adapter decision does not establish what meshers can do. The first table measures microgen’s VTK path at a deliberately small eight points per cell. The second table calls meshers directly on coordinate charts for both volume and surface.</p>
{coverage_table}
<div class="note">All three VTK sample surfaces have open edges at this coarse resolution. Successful generation does not certify them for printing.</div>
<h3>Direct meshers capability probes</h3>
{capability_table}
<div class="note">The full-cylinder volume used 12 points per cell, an identity periodic transform across the angular seam, and a final point weld. It met the adapter’s 0.1 minimum MMG quality and 0.01 sampled-error gates; the welded volume boundary had zero open edges. Full sphere and Sweep volume maps inverted background tetrahedra at their collapsed axes. Both sector volumes passed at 12 points per cell. For surfaces, the mapped cylinder closes only after periodic seam pairing and cleanup. Sphere and Sweep still leave open edges near singular axes.</div>
<p class="muted">These are single-run capability probes, not matched performance comparisons. Surface probes used 24, 48 or 72 uniform cells across each entire chart; VTK used eight points per unit cell. The surface API has no coordinate-map argument, so the probe mapped its generated vertices afterward. The timings omit that mapping and cleanup. The full cylinder needed 72 cells and returned more than 1.3 million triangles, so this is not yet an efficient print-surface replacement.</p>
</section>

<section id="volume"><div class="section-tag">Limit 03 · raw volume speed</div><h2>Raw VTK grids are faster, but serve a different job</h2>
<p>microgen’s legacy volume method returns clipped mixed cells. Meshers returns tetrahedra and measures MMG shape quality and sampled surface error. Raw VTK wins every timing below; these rows should not be read as a matched FEM-ready workflow comparison.</p>
<div class="chart-card">{volume_chart}<div class="legend"><span><i class="swatch" style="background:var(--vtk)"></i>microgen raw VTK</span><span><i class="swatch" style="background:var(--mesh)"></i>meshers tetrahedra</span><span>Right-hand number = meshers / VTK time</span></div></div>
{volume_table}
<div class="note">At 16 grid points, direct split-P meshers volume reached minimum MMG quality {representative("volume", "split_p_unit", "meshers_volume")["minimum_mmg_quality"]:.3f}, below the 0.1 gate used in the microgen adapter. The adapter refines that case. Earlier full microgen + MMG tests favored meshers for the tested FEM cases, but element counts and surface accuracy were not matched, and MMG failed on some larger cases.</div>
<h3>Disabling tetrahedral optimization</h3>
<p>To test the same speed-for-quality tradeoff for volumes, meshers ran with zero optimization passes on three cases. It remains slower than raw VTK’s mixed-cell grid, while its minimum tetrahedral quality falls. Raw VTK cells are not equivalent to FEM-ready tetrahedra, so this is only a lower-bound speed check. The relevant quality-matched volume workflow remains microgen VTK followed by MMG.</p>
{low_quality_volume_table}
</section>

<section id="method"><div class="section-tag">Method and scope</div><h2>How to read the numbers</h2>
<p>The surface speed experiment tests three output paths at similar triangle counts: microgen VTK, meshers without improvement but with exact edge intersections, and an experimental meshers build with linear edge intersections. Its scan considered 10–15 meshers grid points per cell. All timings exclude metric calculations and imports. The linear prototype also uses one field normal per polygon. The fast modes are not the default meshers surface configuration.</p>
<p>The low-optimization volume check used the same 16 grid points per cell for both paths. It ran three fresh-process trials with alternating order. Meshers used <code>optimize_passes=0</code>; geometry tolerance remained 0.01. Cell counts and quality are not matched between raw VTK mixed cells and meshers tetrahedra.</p>
<p>The matched-count runner selects the meshers resolution with the smallest triangle-count gap from 10–15 grid points per cell, against VTK at 16. It then runs three fresh-process trials per path, alternating order. Both paths use the same TPMS shape, thickness and grading. Generation timing stops before triangle metrics. Area and angles use every output triangle; equivalent-area diameter is derived from triangle area. Count matching changes input resolution, so the comparison does not establish equal geometric accuracy.</p>
<p>Each timed path ran in a fresh process, three times, with order alternated. Tables show medians of <code>runtime_seconds</code>, which includes geometry construction and generation but excludes module import and process startup. Both paths used 16 grid points per unit-cell axis. meshers may retry a nearby background resolution for periodic surface quality. Output element counts differ, so these are matched-input-resolution comparisons, not matched-mesh-size or matched-geometric-error comparisons.</p>
<p>Curved VTK examples used eight points per cell to keep the coverage check small. Their runtime includes shape construction. The adapter fallback was checked once per geometry; it has no meshers generation time. The direct surface build was compiled with <code>experimental-surfaces</code>, and microgen’s branch still uses VTK for <code>generate_surface_mesh()</code>.</p>
<p class="muted">Host: {platform_name}. Earlier benchmark base: microgen revision <code>{microgen_rev}</code>; meshers revision <code>{meshers_rev}</code>. Speed-quality probes used experimental meshers revision <code>{prototype_rev}</code>. Source: <a href="fast_surface_benchmark.py">speed-quality runner</a>, <a href="fast_surface_data.json">fast exact JSON</a>, <a href="linear_surface_data.json">linear JSON</a>, <a href="check_fast_surface_topology.py">topology check</a>, <a href="fast_surface_topology.json">topology JSON</a>, <a href="fast_volume_benchmark.py">low-optimization volume runner</a>, <a href="fast_volume_data.json">volume JSON</a>, <a href="matched_surface_benchmark.py">matched-count runner</a>, <a href="matched_surface_data.json">matched-count JSON</a>, <a href="benchmark_limitations.py">fixed-resolution runner</a>, <a href="limitations_data.json">fixed-resolution JSON</a>, <a href="probe_curved_meshers.py">direct curved probe</a>, <a href="curved_core_probe_data.json">curved probe JSON</a>.</p>
<details><summary>Show all surface and volume timing samples</summary><div class="table-wrap"><table><thead><tr><th>Case</th><th>Path</th><th>Trial seconds</th></tr></thead><tbody>{sample_rows}</tbody></table></div></details>
<details><summary>Show curved-geometry timing samples</summary><div class="table-wrap"><table><thead><tr><th>Geometry</th><th>Path</th><th>Trial seconds</th></tr></thead><tbody>{raw_times("curved", CURVED_CASES)}</tbody></table></div></details>
</section>
<footer>Benchmarked on experimental branches. No meshers 0.1.0 surface release or microgen pull request is implied by this report.</footer>
</div></body></html>"""


if __name__ == "__main__":
    output = HERE / "meshers_limitations_report.html"
    output.write_text(build(), encoding="utf-8")
    print(output)
