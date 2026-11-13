# Checks the integral against a simulation of 1M families (Alzheimer's, 2026 inputs, 5 embryos).
import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from model import SITE, expected_risks

K, R2 = SITE["Alzheimer's disease"]
rng = np.random.default_rng(0)
families = 1_000_000
sd_half = np.sqrt(R2 / 2)

scores = rng.normal(0, sd_half, (families, 1)) + rng.normal(0, sd_half, (families, 5))
risk = norm.cdf((scores - norm.isf(K)) / np.sqrt(1 - R2))
best = risk.min(axis=1).mean()
worst = risk.max(axis=1).mean()
avg = risk.mean()

b, a, w = expected_risks(K, R2, 5)
print(f"simulation: vsWorst {1 - best / worst:.2%}  vsAvg {1 - best / avg:.2%}")
print(f"integral:   vsWorst {1 - b / w:.2%}  vsAvg {1 - b / a:.2%}")
