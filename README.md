# Embryo selection risk reductions: an independent check of Nucleus Origin

A small, dependency-free reproduction of the embryo-selection numbers Nucleus Genomics publishes for its nine Origin polygenic scores (Alzheimer's disease, breast cancer, coronary artery disease, endometriosis, hypertension, prostate cancer, rheumatoid arthritis, type 1 diabetes, type 2 diabetes).

**Independent work by Cody Chen. Not affiliated with or endorsed by Nucleus Genomics.** Not medical advice.

## What it does

- `demo/index.html`: open it in a browser (no build step). Pick a disease and 1 to 10 embryos and see the expected risk for the best, average and worst embryo, plus the relative reduction vs worst and vs average. You can switch between the 2026 table and the 2025 preprint inputs, and between within-family and population R². For preprint inputs it also shows a range from the published liability R² 95% CI.
- `demo/calc.js`: the model in about 75 lines of plain JavaScript, no dependencies.
- `check.py`: Python (NumPy/SciPy) version, using different numerics (Gauss-Hermite over the parental midpoint).
- `crosscheck_py.py` + `crosscheck.js`: compares the two implementations over 360 cases (9 diseases x 2 input sets x 2 R² types x 1 to 10 embryos).
- `make_chart.py` -> `chart.png`: vs worst and vs average at 5 embryos for all nine diseases.

## Model

Standard liability threshold model, following Lencz et al. (eLife 2021), which the Nucleus preprint cites:

- Liability = PGS + residual, variance 1. Disease if liability > T = Φ⁻¹(1 − K), where K is lifetime prevalence.
- PGS variance = R² on the liability scale.
- Sibling embryos share the parental midpoint (variance R²/2) and differ by independent deviations (variance R²/2).
- Expected risk of the lowest- and highest-scoring of n embryos comes from integrating over the order-statistic densities. The average embryo's expected risk equals K.

## Results (5 embryos, 2026 table, within-family R²)

| Disease | Prevalence | R² | Best | Average | Worst | vs worst | vs average |
|---|---|---|---|---|---|---|---|
| Alzheimer's disease | 11.1% | 21.1% | 4.96% | 11.1% | 19.2% | 74.2% | 55.3% |
| Breast cancer | 13.7% | 12.9% | 7.85% | 13.7% | 20.8% | 62.2% | 42.7% |
| Coronary artery disease | 40.1% | 21.0% | 25.9% | 40.1% | 55.2% | 53.0% | 35.4% |
| Endometriosis | 11.0% | 6.0% | 7.49% | 11.0% | 15.1% | 50.3% | 31.9% |
| Hypertension | 57.2% | 15.6% | 44.2% | 57.2% | 69.8% | 36.6% | 22.7% |
| Prostate cancer | 12.7% | 21.7% | 5.79% | 12.7% | 21.7% | 73.3% | 54.4% |
| Rheumatoid arthritis | 1.8% | 5.7% | 1.04% | 1.8% | 2.77% | 62.5% | 42.3% |
| Type 1 diabetes | 0.59% | 14.0% | 0.196% | 0.59% | 1.21% | 83.8% | 66.7% |
| Type 2 diabetes | 29.0% | 25.4% | 15.8% | 29.0% | 44.2% | 64.3% | 45.6% |

## Checks

- With the preprint's Supplementary Table 2 inputs, the model reproduces all 27 entries of the preprint's 2, 5 and 10 embryo table to within 1 point.
- With the 2026 table inputs, it matches the expected risks behind the calculator on mynucleus.com to within 0.05 points (rounding of prevalence).
- JS and Python agree to within 1e-10 percentage points across all 360 cases.
- Monte Carlo (1M families) agrees with the integral to 0.05 points.

## Observation

The preprint's vs-average range (27 to 67%, including 50% for type 2 diabetes) uses the preprint's prevalence assumptions (type 2 diabetes 19.9%, hypertension 46.85%). The 2026 table uses 29% and 57.2%, and its R² values match the same models rescaled to those prevalences. With the 2026 inputs, type 2 diabetes becomes 45.6% and hypertension 22.7%, so the five-embryo range is about 23 to 67%. The other seven diseases move by less than 1 point.

## Sources

- Nucleus Labs / Genetic Optimization Hub, "UK Biobank Validation of Nucleus Models": https://mynucleus.com/labs (read Nov 12, 2026)
- Nucleus embryo page: https://mynucleus.com/deep-dives/embryo
- Cordogan et al., medRxiv 10.1101/2025.10.24.25338613 (v1 to v3; Supplementary Tables 2 and 3 are the same in all versions)
- Lencz T. et al., "Utility of polygenic embryo screening for disease depends on the selection strategy", eLife 2021

## Run

```
python3 check.py              # tables + sensitivity + Monte Carlo
python3 crosscheck_py.py && node crosscheck.js
python3 make_chart.py
```
