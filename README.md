# Checking Nucleus Origin's embryo-selection numbers

I was reading Nucleus Genomics' preprint on their Origin polygenic scores and wanted to see if I could rebuild their embryo-selection numbers myself. The model reproduces their published figures, but I found that the embryo page's type 2 diabetes and hypertension figures use older prevalences than their 2026 UK Biobank table. This is my own side project and I'm not affiliated with Nucleus. It's not medical advice.

## What's here

- `demo/index.html` is a small calculator you can open in a browser. Pick a disease and 1 to 10 embryos to see the expected risk for the best, average and worst embryo. You can switch between the 2026 table and the 2025 preprint inputs.
- `demo/calc.js` is the model in JavaScript.
- `model.py` is the same model in Python (NumPy/SciPy), plus the inputs from both sources.
- `reproduce.py` prints the match against the preprint's table and the 2026 type 2 diabetes and hypertension comparison.
- `tests/` compares the JavaScript and Python versions and checks the integral against a simulation.
- `scripts/make_chart.py` draws `chart.png`.

## The model

I used the liability threshold model from Lencz et al. (eLife 2021), which the Nucleus preprint cites. Liability is the polygenic score plus a residual, with total variance 1, and you get the disease if liability is above the threshold set by lifetime prevalence. The score's variance is R² on the liability scale. Sibling embryos share the parents' midpoint (variance R²/2) and differ by independent deviations (variance R²/2). The best and worst of n embryos come from integrating over order statistics, and the average embryo's risk is just the prevalence.

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
| Type 2 diabetes | 29.0% | 25.4% | 15.8% | 29.0% | 44.2% | 64.3% | 45.5% |

The type 2 diabetes vs average figure is 45.5%, which is the value Nucleus's calculator shows. My model gives 45.6%, and the difference comes from rounding of the prevalence.

## What I found

The preprint gives a vs average range of 27 to 67%, including 50% for type 2 diabetes. Those numbers use the preprint's prevalences (19.9% for type 2 diabetes and 46.85% for hypertension). The 2026 table uses 29% and 57.2%, and its R² values match the same models rescaled to those prevalences. With the 2026 inputs, type 2 diabetes comes out at 45.5% and hypertension at 22.7%, so the five-embryo range is about 23 to 67%. The other seven diseases move by less than 1 point.

## How I checked it

With the preprint's Supplementary Table 2 inputs, the model reproduces all 27 entries of the preprint's 2, 5 and 10 embryo table to within 1 point. With the 2026 inputs, it matches the calculator on mynucleus.com to within 0.05 points. The JavaScript and Python versions agree to within 0.001 points across 180 cases, and a Monte Carlo run with 1M families agrees with the integral to 0.05 points.

- Nucleus validation table: https://mynucleus.com/labs (read Nov 12, 2026)
- Nucleus embryo page: https://mynucleus.com/deep-dives/embryo
- Cordogan et al., medRxiv 10.1101/2025.10.24.25338613

## Running it

```
python3 reproduce.py
python3 tests/crosscheck_py.py && node tests/crosscheck.js
python3 tests/monte_carlo.py
python3 scripts/make_chart.py
```
