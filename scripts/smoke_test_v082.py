from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


manifest_errors = []
with (PACKAGE / "MANIFEST_SHA256.csv").open(encoding="utf-8-sig", newline="") as stream:
    manifest = list(csv.DictReader(stream))
for row in manifest:
    path = PACKAGE / row["relative_path"]
    if not path.exists() or path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
        manifest_errors.append(row["relative_path"])
if manifest_errors:
    raise SystemExit(json.dumps({"status": "FAIL_MANIFEST", "errors": manifest_errors}))

generated = PACKAGE / "generated"
if generated.exists():
    shutil.rmtree(generated)
for name in ("build_publication_figures_v082.py", "build_supplementary_figures_v082.py"):
    subprocess.run([sys.executable, str(PACKAGE / "scripts" / name)], cwd=PACKAGE, check=True, timeout=300)

errors = []
with (PACKAGE / "EXPECTED_FIGURE_HASHES.csv").open(encoding="utf-8-sig", newline="") as stream:
    expected = list(csv.DictReader(stream))
for row in expected:
    path = PACKAGE / row["relative_path"]
    if not path.exists() or path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
        errors.append(row["relative_path"])
result = {
    "status": "PASS" if not errors else "FAIL_OUTPUT_HASH",
    "manifest_entries": len(manifest),
    "expected_outputs": len(expected),
    "mismatches": errors,
}
print(json.dumps(result))
raise SystemExit(0 if not errors else 1)