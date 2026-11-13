// Compares demo/calc.js with the Python reference values (difference in percentage points).
const calc = require("../demo/calc.js");
const ref = require("./crosscheck_ref.json");

let maxDiff = 0;
let worstCase = "";
for (const row of ref) {
  const disease = calc.DISEASES.find((d) => d.name === row.name);
  const inp = calc.inputs(disease, row.src);
  const js = calc.expectedRisks(inp.K, inp.R2, row.n);
  for (const key of ["best", "avg", "worst", "vsWorst", "vsAvg"]) {
    const diff = Math.abs(100 * (js[key] - row[key]));
    if (diff > maxDiff) {
      maxDiff = diff;
      worstCase = `${row.name} ${row.src} n=${row.n} ${key}`;
    }
  }
}
console.log(`rows ${ref.length}, max diff ${maxDiff.toExponential(2)} points (${worstCase})`);

console.log("2026 table, 5 embryos:");
for (const d of calc.DISEASES) {
  const inp = calc.inputs(d, "site");
  const r = calc.expectedRisks(inp.K, inp.R2, 5);
  const pct = (x) => (100 * x).toFixed(2) + "%";
  console.log(`  ${d.name.padEnd(24)} best ${pct(r.best)}  avg ${pct(r.avg)}  worst ${pct(r.worst)}  vsWorst ${(100 * r.vsWorst).toFixed(1)}%  vsAvg ${(100 * r.vsAvg).toFixed(1)}%`);
}

process.exit(maxDiff < 0.1 ? 0 : 1);
