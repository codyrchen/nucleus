import numpy as np
from scipy.stats import norm

# Grid for the parents' midpoint and for each embryo's own deviation
GH_X, GH_W = np.polynomial.hermite_e.hermegauss(80)
GH_W = GH_W / GH_W.sum()
Z = np.linspace(-9, 9, 3601)
DZ = Z[1] - Z[0]


def expected_risks(K, R2, n=5):
    """Expected risk of the best, average and worst of n embryos."""
    T = norm.isf(K)
    sd_rest = np.sqrt(1 - R2)
    sd_half = np.sqrt(R2 / 2)
    pdf = norm.pdf(Z)
    f_best = n * pdf * norm.sf(Z) ** (n - 1)
    f_worst = n * pdf * norm.cdf(Z) ** (n - 1)
    m = sd_half * GH_X[:, None]
    risk = norm.cdf((m + sd_half * Z - T) / sd_rest)
    best = GH_W @ (risk @ f_best) * DZ
    avg = GH_W @ (risk @ pdf) * DZ
    worst = GH_W @ (risk @ f_worst) * DZ
    return best, avg, worst


# 2025 preprint: (prevalence, R2). T1D R2 is scaled down for within-family selection.
PAPER = {
    "Alzheimer's disease": (0.111, 0.211),
    "Breast cancer": (0.1387, 0.128),
    "Coronary artery disease": (0.405, 0.210),
    "Endometriosis": (0.110, 0.060),
    "Hypertension": (0.4685, 0.155),
    "Prostate cancer": (0.1257, 0.216),
    "Rheumatoid arthritis": (0.026, 0.064),
    "Type 1 diabetes": (0.0067, 0.210 * 0.827**2),
    "Type 2 diabetes": (0.199, 0.229),
}

# 2026 table: (prevalence, within-family R2)
SITE = {
    "Alzheimer's disease": (0.111, 0.211),
    "Breast cancer": (0.137, 0.129),
    "Coronary artery disease": (0.401, 0.210),
    "Endometriosis": (0.110, 0.060),
    "Hypertension": (0.572, 0.156),
    "Prostate cancer": (0.127, 0.217),
    "Rheumatoid arthritis": (0.018, 0.057),
    "Type 1 diabetes": (0.0059, 0.140),
    "Type 2 diabetes": (0.290, 0.254),
}
