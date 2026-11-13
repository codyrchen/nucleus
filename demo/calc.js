// Liability threshold embryo-selection calculator. Independent reproduction by Cody Chen.
// Not affiliated with Nucleus Genomics.
(function (root) {
  // erfc: Taylor series for erf when |x| < 2.5, backward continued fraction otherwise (~1e-13 abs error)
  function erfc(x) {
    const z = Math.abs(x);
    let r;
    if (z < 2.5) {
      let sum = z, term = z, n = 0;
      do { n++; term *= -z * z / n; sum += term / (2 * n + 1); } while (Math.abs(term) > 1e-17 * Math.abs(sum));
      r = 1 - 2 / Math.sqrt(Math.PI) * sum;
    } else {
      // erfc(z) = exp(-z^2)/sqrt(pi) / (z + (1/2)/(z + (2/2)/(z + (3/2)/(z + ...))))
      let f = z;
      for (let k = 120; k >= 1; k--) f = z + (k / 2) / f;
      r = Math.exp(-z * z) / Math.sqrt(Math.PI) / f;
    }
    return x >= 0 ? r : 2 - r;
  }
  const Phi = (x) => 0.5 * erfc(-x / Math.SQRT2);
  const phi = (x) => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
  function invPhi(p) { // Acklam + two Newton steps
    const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239],
      b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572],
      c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783],
      d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
    let x;
    if (p < 0.02425) { const q = Math.sqrt(-2 * Math.log(p)); x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
    else if (p <= 0.97575) { const q = p - 0.5, r = q*q; x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1); }
    else { const q = Math.sqrt(-2 * Math.log(1 - p)); x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
    for (let i = 0; i < 2; i++) x -= (Phi(x) - p) / phi(x);
    return x;
  }
  // Expected risk of best (lowest-PGS), random, and worst (highest-PGS) of n sibling embryos.
  // Embryo PGS = m + d, m,d ~ N(0, R2/2). Ranking depends only on d; integrating m analytically gives
  // P(disease | d) = Phi((sqrt(R2/2) d - T) / sqrt(1 - R2/2)). Then integrate over the order-statistic density of d.
  function expectedRisks(K, R2, n) {
    const T = invPhi(1 - K), sw = Math.sqrt(R2 / 2), se = Math.sqrt(1 - R2 / 2);
    const N = 3601, lo = -9, hi = 9, h = (hi - lo) / (N - 1);
    let best = 0, worst = 0, avg = 0;
    for (let i = 0; i < N; i++) {
      const z = lo + i * h, w = (i === 0 || i === N - 1) ? 0.5 : 1;
      const p = phi(z), F = Phi(z), r = Phi((sw * z - T) / se);
      best += w * n * p * Math.pow(1 - F, n - 1) * r;
      worst += w * n * p * Math.pow(F, n - 1) * r;
      avg += w * p * r;
    }
    best *= h; worst *= h; avg *= h;
    return { best, avg, worst, vsWorst: 1 - best / worst, vsAvg: 1 - best / avg };
  }
  // Data. SITE: "UK Biobank Validation of Nucleus Models", mynucleus.com/labs (read 2026-11-12).
  // PAPER: medRxiv 10.1101/2025.10.24.25338613 (v1-v3 identical), Supp. Table 2 (liability R2 + 95% CI)
  // and Supp. Table 3 (within-family / population beta ratio).
  const DISEASES = [
    // name, site K, site popR2, site wfR2, paper K, paper R2, paper CI lo, hi, beta ratio
    ["Alzheimer's disease",     0.111,  0.211, 0.211, 0.111,  0.211, 0.148, 0.274, 0.873],
    ["Breast cancer",           0.137,  0.129, 0.129, 0.1387, 0.128, 0.111, 0.146, 0.911],
    ["Coronary artery disease", 0.401,  0.210, 0.210, 0.405,  0.210, 0.186, 0.233, 1.014],
    ["Endometriosis",           0.110,  0.060, 0.060, 0.110,  0.060, 0.042, 0.078, 1.026],
    ["Hypertension",            0.572,  0.156, 0.156, 0.4685, 0.155, 0.146, 0.163, 0.978],
    ["Prostate cancer",         0.127,  0.217, 0.217, 0.1257, 0.216, 0.187, 0.245, 1.025],
    ["Rheumatoid arthritis",    0.018,  0.057, 0.057, 0.026,  0.064, 0.044, 0.083, 0.967],
    ["Type 1 diabetes",         0.0059, 0.204, 0.140, 0.0067, 0.210, 0.144, 0.275, 0.827],
    ["Type 2 diabetes",         0.290,  0.254, 0.254, 0.199,  0.229, 0.211, 0.247, 0.994],
  ].map(([name, sK, sP, sW, pK, pR, pLo, pHi, ratio]) => ({ name, sK, sP, sW, pK, pR, pLo, pHi, ratio }));
  // Returns {K, R2, ci:[lo,hi]|null} for a given source and R2 type.
  // Paper convention (Results section): only T1D within-family differs significantly, so only T1D R2 is
  // reduced (by beta ratio squared) for within-family selection.
  function inputs(d, source, wf) {
    if (source === "site") return { K: d.sK, R2: wf ? d.sW : d.sP, ci: null };
    const adj = (wf && d.name === "Type 1 diabetes") ? d.ratio * d.ratio : 1;
    return { K: d.pK, R2: d.pR * adj, ci: [d.pLo * adj, d.pHi * adj] };
  }
  const api = { erfc, Phi, invPhi, expectedRisks, DISEASES, inputs };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.EmbryoCalc = api;
})(this);
