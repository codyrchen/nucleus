from scipy.stats import norm

from model import PAPER, SITE, expected_risks

# Published "vs average" reductions at 5 embryos from the preprint
PUBLISHED_5 = {
    "Alzheimer's disease": 0.55,
    "Breast cancer": 0.42,
    "Coronary artery disease": 0.35,
    "Endometriosis": 0.31,
    "Hypertension": 0.27,
    "Prostate cancer": 0.54,
    "Rheumatoid arthritis": 0.42,
    "Type 1 diabetes": 0.67,
    "Type 2 diabetes": 0.50,
}


def reductions(K, R2, n=5):
    """Relative risk reduction of the best embryo vs worst and vs average."""
    best, avg, worst = expected_risks(K, R2, n)
    return 1 - best / worst, 1 - best / avg


def liability_scale(K):
    """How liability R2 scales with prevalence for a fixed observed-scale R2."""
    return K**2 * (1 - K) ** 2 / norm.pdf(norm.isf(K)) ** 2


def main():
    print("Preprint inputs vs the preprint's table (vs average)")
    print(f"{'disease':25s} {'n=2':>6s} {'n=5':>6s} {'pub':>5s} {'n=10':>6s}")
    for name, (K, R2) in PAPER.items():
        n2 = reductions(K, R2, 2)[1]
        n5 = reductions(K, R2, 5)[1]
        n10 = reductions(K, R2, 10)[1]
        print(f"{name:25s} {n2:6.1%} {n5:6.1%} {PUBLISHED_5[name]:5.0%} {n10:6.1%}")

    print("\nType 2 diabetes and hypertension, 5 embryos, vs average")
    for name in ["Type 2 diabetes", "Hypertension"]:
        K_old, R2_old = PAPER[name]
        K_new, R2_new = SITE[name]
        print(f"  {name}")
        print(f"    preprint inputs:      {reductions(K_old, R2_old)[1]:.1%}")
        print(f"    2026 prevalence only: {reductions(K_new, R2_old)[1]:.1%}")
        print(f"    2026 R2 only:         {reductions(K_old, R2_new)[1]:.1%}")
        print(f"    2026 inputs:          {reductions(K_new, R2_new)[1]:.1%}")
        rescaled = R2_old * liability_scale(K_new) / liability_scale(K_old)
        print(f"    preprint R2 rescaled to 2026 prevalence: {rescaled:.3f} (2026 table: {R2_new:.3f})")


if __name__ == "__main__":
    main()
