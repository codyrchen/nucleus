# Draws chart.png: best of 5 embryos vs worst and vs average, 2026 inputs.
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from model import SITE
from reproduce import PUBLISHED_5, reductions

names = list(SITE)
short = ["Alzheimer's", "Breast ca.", "CAD", "Endometriosis", "Hypertension", "Prostate ca.", "RA", "T1D", "T2D"]
vs_worst = np.array([reductions(*SITE[d], 5)[0] for d in names]) * 100
vs_avg = np.array([reductions(*SITE[d], 5)[1] for d in names]) * 100
published = np.array([PUBLISHED_5[d] for d in names]) * 100

x = np.arange(len(names))
w = 0.38
fig, ax = plt.subplots(figsize=(10, 5.2), dpi=150)
bars_worst = ax.bar(x - w / 2, vs_worst, w, label="vs worst embryo (what the calculator shows)", color="#d0533f")
bars_avg = ax.bar(x + w / 2, vs_avg, w, label="vs average embryo", color="#2f6fdf")
ax.scatter(x + w / 2, published, marker="_", s=420, color="black", linewidths=2, zorder=3, label="Published vs average (2025 preprint)")
for bar, v in zip(bars_worst, vs_worst):
    ax.text(bar.get_x() + bar.get_width() / 2, v + 1, f"{v:.0f}", ha="center", fontsize=8)
for bar, v in zip(bars_avg, vs_avg):
    ax.text(bar.get_x() + bar.get_width() / 2, v + 1, f"{v:.0f}", ha="center", fontsize=8)

ax.set_xticks(x, short, fontsize=9)
ax.set_ylabel("Relative risk reduction, best of 5 embryos (%)")
ax.set_ylim(0, 95)
ax.set_title("Nucleus Origin, 5 embryos: lowest-risk embryo vs worst and vs average\n(2026 UK Biobank table inputs, within-family R², liability threshold model)", fontsize=11)
ax.legend(fontsize=8, loc="upper left", frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.text(0.01, 0.01, "Independent reproduction by Cody Chen, not affiliated with Nucleus Genomics. Sources: mynucleus.com/labs; medRxiv 10.1101/2025.10.24.25338613.", fontsize=7, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(ROOT / "chart.png")
print("saved chart.png")
