import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
from check import rrr, SITE, PAPER_PUBLISHED_5
names = list(SITE); short = ["Alzheimer's","Breast ca.","CAD","Endometriosis","Hypertension","Prostate ca.","RA","T1D","T2D"]
vw, va = zip(*[rrr(SITE[d][0], SITE[d][2], 5)[:2] for d in names])
pub = [PAPER_PUBLISHED_5[d] for d in names]
x = np.arange(len(names)); w = 0.38
fig, ax = plt.subplots(figsize=(10, 5.2), dpi=150)
b1 = ax.bar(x - w/2, np.array(vw)*100, w, label="vs worst embryo (what the calculator shows)", color="#d0533f")
b2 = ax.bar(x + w/2, np.array(va)*100, w, label="vs average embryo", color="#2f6fdf")
ax.scatter(x + w/2, np.array(pub)*100, marker="_", s=420, color="black", linewidths=2, zorder=3, label="Published vs average (2025 preprint)")
for r, v in zip(b1, vw): ax.text(r.get_x()+r.get_width()/2, v*100+1, f"{v*100:.0f}", ha="center", fontsize=8)
for r, v in zip(b2, va): ax.text(r.get_x()+r.get_width()/2, v*100+1, f"{v*100:.0f}", ha="center", fontsize=8)
ax.set_xticks(x, short, fontsize=9); ax.set_ylabel("Relative risk reduction, best of 5 embryos (%)"); ax.set_ylim(0, 95)
ax.set_title("Nucleus Origin, 5 embryos: lowest-risk embryo vs worst and vs average\n(2026 UK Biobank table inputs, within-family R², liability threshold model)", fontsize=11)
ax.legend(fontsize=8, loc="upper left", frameon=False); ax.spines[["top","right"]].set_visible(False)
fig.text(0.01, 0.01, "Independent reproduction by Cody Chen, not affiliated with Nucleus Genomics. Sources: mynucleus.com/labs; medRxiv 10.1101/2025.10.24.25338613.", fontsize=7, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig("chart.png")
print("saved chart.png")
