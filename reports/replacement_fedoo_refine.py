import json
import subprocess
import sys
from pathlib import Path

p = Path("reports/replacement_fedoo.json")
rows = json.loads(p.read_text())
old = list(rows)
for mode, n in [("mmg", 16), ("mmg", 20), ("meshers", 36)]:
    r = subprocess.run(
        [sys.executable, "reports/replacement_fedoo.py", mode, str(n)],
        capture_output=True,
        text=True,
        timeout=180,
    )
    row = (
        dict(json.loads(r.stdout.strip().splitlines()[-1]), status="ok")
        if r.returncode == 0
        else dict(mode=mode, resolution=n, status="error", error=r.stderr[-2500:])
    )
    rows = [x for x in rows if (x["mode"], x["resolution"]) != (mode, n)] + [row]
    p.write_text(json.dumps(rows, indent=2))
    print(json.dumps(row), flush=True)
Path("reports/replacement_fedoo_initial.json").write_text(json.dumps(old, indent=2))
