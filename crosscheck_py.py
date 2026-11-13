# Writes Python reference values (check.py model) for the JS cross-check.
import json
from check import expected_risks, PAPER, SITE
out = []
for d in PAPER:
    for src in ("site", "paper"):
        for wf in (True, False):
            if src == "site": K = SITE[d][0]; R2 = SITE[d][2] if wf else SITE[d][1]
            else:
                K, R2, ci, ratio = PAPER[d]
                if wf and d == "Type 1 diabetes": R2 *= ratio ** 2
            for n in range(1, 11):
                b, a, w = expected_risks(K, R2, n)
                out.append(dict(name=d, src=src, wf=wf, n=n, best=b, avg=a, worst=w,
                                vsWorst=1 - b / w, vsAvg=1 - b / a))
json.dump(out, open("crosscheck_ref.json", "w"))
print(len(out), "reference rows")
