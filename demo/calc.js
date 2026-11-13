// Embryo selection calculator using a liability threshold model.
(function (root) {
  // Standard normal density
  function phi(x) {
    return Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
  }

  // Standard normal CDF (polynomial approximation, error below 1e-7)
  function Phi(x) {
    const t = 1 / (1 + 0.2316419 * Math.abs(x));
    const poly = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    const upper = phi(x) * poly;
    return x >= 0 ? 1 - upper : upper;
  }

  // Inverse of Phi by bisection
  function invPhi(p) {
    let lo = -10;
    let hi = 10;
    for (let i = 0; i < 100; i++) {
      const mid = (lo + hi) / 2;
      if (Phi(mid) < p) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    return (lo + hi) / 2;
  }

  // Expected risk of the best, average and worst of n sibling embryos.
  // Each embryo's score is a shared parental part plus its own deviation, each with variance R2/2.
  // Ranking only depends on the deviation, so the shared part can be folded into the risk.
  function expectedRisks(K, R2, n) {
    const T = invPhi(1 - K);
    const sdHalf = Math.sqrt(R2 / 2);
    const sdRest = Math.sqrt(1 - R2 / 2);
    const steps = 3601;
    const lo = -9;
    const h = 18 / (steps - 1);
    let best = 0;
    let worst = 0;
    let avg = 0;
    for (let i = 0; i < steps; i++) {
      const z = lo + i * h;
      const weight = i === 0 || i === steps - 1 ? 0.5 : 1;
      const p = phi(z);
      const F = Phi(z);
      const risk = Phi((sdHalf * z - T) / sdRest);
      best += weight * n * p * Math.pow(1 - F, n - 1) * risk;
      worst += weight * n * p * Math.pow(F, n - 1) * risk;
      avg += weight * p * risk;
    }
    best *= h;
    worst *= h;
    avg *= h;
    return { best, avg, worst, vsWorst: 1 - best / worst, vsAvg: 1 - best / avg };
  }

  // name, 2026 prevalence, 2026 within-family R2, preprint prevalence, preprint R2
  // T1D preprint R2 is scaled down for within-family selection.
  const DISEASES = [
    ["Alzheimer's disease", 0.111, 0.211, 0.111, 0.211],
    ["Breast cancer", 0.137, 0.129, 0.1387, 0.128],
    ["Coronary artery disease", 0.401, 0.210, 0.405, 0.210],
    ["Endometriosis", 0.110, 0.060, 0.110, 0.060],
    ["Hypertension", 0.572, 0.156, 0.4685, 0.155],
    ["Prostate cancer", 0.127, 0.217, 0.1257, 0.216],
    ["Rheumatoid arthritis", 0.018, 0.057, 0.026, 0.064],
    ["Type 1 diabetes", 0.0059, 0.140, 0.0067, 0.210 * 0.827 * 0.827],
    ["Type 2 diabetes", 0.290, 0.254, 0.199, 0.229],
  ].map(([name, siteK, siteR2, paperK, paperR2]) => ({ name, siteK, siteR2, paperK, paperR2 }));

  // Prevalence and R2 for the chosen input set ("site" or "paper")
  function inputs(d, source) {
    if (source === "site") {
      return { K: d.siteK, R2: d.siteR2 };
    }
    return { K: d.paperK, R2: d.paperR2 };
  }

  const api = { Phi, invPhi, expectedRisks, DISEASES, inputs };
  if (typeof module !== "undefined" && module.exports) {
    module.exports = api;
  } else {
    root.EmbryoCalc = api;
  }
})(this);
