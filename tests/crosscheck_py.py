# Writes Python reference values for the JavaScript cross-check.
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from model import PAPER, SITE, expected_risks

rows = []
for source, table in [("site", SITE), ("paper", PAPER)]:
    for name, (K, R2) in table.items():
        for n in range(1, 11):
            best, avg, worst = expected_risks(K, R2, n)
            rows.append({
                "name": name,
                "src": source,
                "n": n,
                "best": best,
                "avg": avg,
                "worst": worst,
                "vsWorst": 1 - best / worst,
                "vsAvg": 1 - best / avg,
            })

out = Path(__file__).parent / "crosscheck_ref.json"
out.write_text(json.dumps(rows))
print(len(rows), "reference rows")
