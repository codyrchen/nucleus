"""Independent reproduction of Nucleus Origin embryo-selection risk reductions.
Model (Lencz et al. 2021 style liability threshold):
  liability L = G + E ~ N(0,1); PGS G ~ N(0, R2) with R2 on the liability scale.
  Sibling embryos: G_i = m + d_i, m ~ N(0, R2/2) (parental midpoint), d_i ~ N(0, R2/2) iid.
  Disease if L > T, T = Phi^-1(1-K).  P(disease | G) = Phi((G - T)/sqrt(1 - R2)).
Outputs, for n embryos: expected risk of lowest-PGS embryo (best), random embryo (= K),
highest-PGS embryo (worst); RRR vs worst = 1 - best/worst, RRR vs average = 1 - best/K.
Numerics: Gauss-Hermite over m, fine trapezoid grid over the order-statistic density of d.
Author: Cody Chen. Not affiliated with Nucleus Genomics.
"""
import numpy as np
from scipy.stats import norm

GH_X, GH_W = np.polynomial.hermite_e.hermegauss(80)   # probabilists' Hermite
GH_W = GH_W / GH_W.sum()
Z = np.linspace(-9, 9, 3601); DZ = Z[1] - Z[0]

def expected_risks(K, R2, n=5):
    T = norm.isf(K); se = np.sqrt(1 - R2); sw = np.sqrt(R2 / 2)
    pdf, cdf, sf = norm.pdf(Z), norm.cdf(Z), norm.sf(Z)
    f_min = n * pdf * sf ** (n - 1); f_max = n * pdf * cdf ** (n - 1)
    m = sw * GH_X[:, None]
    r = norm.cdf((m + sw * Z[None, :] - T) / se)          # risk grid [m, z]
    best = GH_W @ (r @ f_min) * DZ
    worst = GH_W @ (r @ f_max) * DZ
    avg = GH_W @ (r @ pdf) * DZ
    return best, avg, worst

def rrr(K, R2, n=5):
    b, a, w = expected_risks(K, R2, n)
    return 1 - b / w, 1 - b / a, b, a, w

# Source A: medRxiv 10.1101/2025.10.24.25338613, Supp. Table 2 (identical in v1/v2/v3)
#   (lifetime prevalence, liability R2, liability R2 95% CI, within-family beta ratio from Supp. Table 3)
PAPER = {
 "Alzheimer's disease":     (0.111,  0.211, (0.148, 0.274), 0.873),
 "Breast cancer":           (0.1387, 0.128, (0.111, 0.146), 0.911),
 "Coronary artery disease": (0.405,  0.210, (0.186, 0.233), 1.014),
 "Endometriosis":           (0.110,  0.060, (0.042, 0.078), 1.026),
 "Hypertension":            (0.4685, 0.155, (0.146, 0.163), 0.978),
 "Prostate cancer":         (0.1257, 0.216, (0.187, 0.245), 1.025),
 "Rheumatoid arthritis":    (0.026,  0.064, (0.044, 0.083), 0.967),
 "Type 1 diabetes":         (0.0067, 0.210, (0.144, 0.275), 0.827),
 "Type 2 diabetes":         (0.199,  0.229, (0.211, 0.247), 0.994),
}
PAPER_PUBLISHED_5 = {"Alzheimer's disease": .55, "Breast cancer": .42, "Coronary artery disease": .35,
  "Endometriosis": .31, "Hypertension": .27, "Prostate cancer": .54, "Rheumatoid arthritis": .42,
  "Type 1 diabetes": .67, "Type 2 diabetes": .50}
# Source B: "UK Biobank Validation of Nucleus Models", mynucleus.com/labs (fetched 2026-11-12)
#   (lifetime prevalence, population R2, within-family R2)
SITE = {
 "Alzheimer's disease":     (0.111,  0.211, 0.211),
 "Breast cancer":           (0.137,  0.129, 0.129),
 "Coronary artery disease": (0.401,  0.210, 0.210),
 "Endometriosis":           (0.110,  0.060, 0.060),
 "Hypertension":            (0.572,  0.156, 0.156),
 "Prostate cancer":         (0.127,  0.217, 0.217),
 "Rheumatoid arthritis":    (0.018,  0.057, 0.057),
 "Type 1 diabetes":         (0.0059, 0.204, 0.140),
 "Type 2 diabetes":         (0.290,  0.254, 0.254),
}

def paper_r2_for_selection(name):
    K, R2, ci, ratio = PAPER[name]
    # paper: only T1D within-family differs significantly; "Liability R2 ... reduced accordingly"
    return R2 * ratio ** 2 if name == "Type 1 diabetes" else R2

if __name__ == "__main__":
    print("A) Paper inputs (Supp. Table 2), 5 embryos, vs average")
    print(f"{'disease':25s} {'K':>7s} {'R2':>6s} | {'mine':>6s} {'pub':>5s} | {'vsWorst':>7s}")
    for d in PAPER:
        K = PAPER[d][0]; R2 = paper_r2_for_selection(d)
        vw, va, *_ = rrr(K, R2, 5)
        print(f"{d:25s} {K:7.2%} {R2:6.3f} | {va:6.1%} {PAPER_PUBLISHED_5[d]:5.0%} | {vw:7.1%}")
    print("\n  2 and 10 embryos (paper inputs) vs average:")
    for d in PAPER:
        K = PAPER[d][0]; R2 = paper_r2_for_selection(d)
        print(f"  {d:25s} n=2 {rrr(K,R2,2)[1]:5.1%}  n=10 {rrr(K,R2,10)[1]:5.1%}")
    print("\nB) 2026 site table inputs, 5 embryos (within-family R2)")
    print(f"{'disease':25s} {'K':>7s} {'R2wf':>6s} | {'best':>6s} {'avg':>6s} {'worst':>6s} | {'vsWorst':>7s} {'vsAvg':>6s}")
    for d, (K, Rp, Rw) in SITE.items():
        vw, va, b, a, w = rrr(K, Rw, 5)
        print(f"{d:25s} {K:7.2%} {Rw:6.3f} | {b:6.2%} {a:6.2%} {w:6.2%} | {vw:7.1%} {va:6.1%}")
    print("\nC) Isolate T2D / hypertension drivers (vs average, 5 embryos)")
    for d in ["Type 2 diabetes", "Hypertension"]:
        Kp, Rp = PAPER[d][0], PAPER[d][1]; Ks, Rs = SITE[d][0], SITE[d][2]
        print(f"  {d}: paper K,R2 {rrr(Kp,Rp)[1]:.1%} | site K only {rrr(Ks,Rp)[1]:.1%} | "
              f"site R2 only {rrr(Kp,Rs)[1]:.1%} | site both {rrr(Ks,Rs)[1]:.1%}")
    # Does the site R2 equal the paper's model re-scaled to the new prevalence?
    # Lee et al. 2012: liability R2 scales with K^2(1-K)^2 / phi(t)^2 for fixed observed-scale R2
    # (ascertainment term ignored; approximate).
    f = lambda K: K**2 * (1-K)**2 / norm.pdf(norm.isf(K))**2
    for d in ["Type 2 diabetes", "Hypertension", "Rheumatoid arthritis", "Type 1 diabetes"]:
        Kp, Rp = PAPER[d][0], PAPER[d][1]; Ks, Rs = SITE[d][0], SITE[d][1]
        print(f"  {d}: paper R2 rescaled to site K ~ {Rp*f(Ks)/f(Kp):.3f} (site lists {Rs:.3f})")
    print("\nD) Alzheimer's sensitivity to the paper's R2 95% CI (14.8-27.4%), 5 embryos")
    for r in (0.148, 0.211, 0.274):
        vw, va, *_ = rrr(0.111, r, 5); print(f"  R2={r:.3f}: vsWorst {vw:.1%}  vsAvg {va:.1%}")
    # Monte Carlo sanity check (AD, site inputs)
    rng = np.random.default_rng(0); N = 1_000_000; K, R = 0.111, 0.211; sw = np.sqrt(R/2)
    G = rng.normal(0, sw, (N, 1)) + rng.normal(0, sw, (N, 5))
    risk = norm.cdf((G - norm.isf(K)) / np.sqrt(1 - R))
    print(f"\nMC check AD: vsWorst {1-risk.min(1).mean()/risk.max(1).mean():.2%} vsAvg {1-risk.min(1).mean()/risk.mean():.2%}")
