// Compares demo/calc.js against Python reference values (max abs diff in percentage points).
const c = require("./demo/calc.js");
const ref = require("./crosscheck_ref.json");
let maxd = 0, worstRow = null;
for (const r of ref) {
  const d = c.DISEASES.find((x) => x.name === r.name);
  const i = c.inputs(d, r.src, r.wf);
  const j = c.expectedRisks(i.K, i.R2, r.n);
  for (const k of ["best", "avg", "worst", "vsWorst", "vsAvg"]) {
    const diff = Math.abs(100 * (j[k] - r[k]));
    if (diff > maxd) { maxd = diff; worstRow = { ...r, key: k, js: j[k] }; }
  }
}
console.log(`rows ${ref.length}, max abs diff ${maxd.toExponential(2)} pt`, worstRow && `${worstRow.name} ${worstRow.src} n=${worstRow.n} ${worstRow.key}`);
console.log("site, within-family, n=5:");
for (const d of c.DISEASES) {
  const i = c.inputs(d, "site", true), j = c.expectedRisks(i.K, i.R2, 5);
  console.log(`  ${d.name.padEnd(24)} best ${(100*j.best).toFixed(2)}%  avg ${(100*j.avg).toFixed(2)}%  worst ${(100*j.worst).toFixed(2)}%  vsWorst ${(100*j.vsWorst).toFixed(1)}%  vsAvg ${(100*j.vsAvg).toFixed(1)}%`);
}
process.exit(maxd < 0.1 ? 0 : 1);
