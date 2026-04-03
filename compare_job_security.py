"""
compare_job_security.py

Compares two employment security models:
  - US-style: at-will employment with sudden termination
  - Japan-style: lifetime employment culture with ongoing termination threat

Each model is evaluated across several dimensions. A score from 0.0 (worst) to
1.0 (best) is assigned per dimension, then an overall weighted score is computed
to determine which model performs better overall.
"""

DIMENSIONS = {
    "employee_stability":   "Stability of employment (predictable tenure)",
    "company_flexibility":  "Employer ability to restructure quickly",
    "chronic_stress":       "Low chronic stress for employees (higher = less stress)",
    "market_dynamism":      "Active, fluid labor market enabling career mobility",
    "career_clarity":       "Employee clarity on current standing / performance",
    "economic_efficiency":  "Efficient reallocation of labor to productive roles",
}

# Equal weights across all dimensions
WEIGHTS = {dim: 1.0 / len(DIMENSIONS) for dim in DIMENSIONS}


def score_us_employment():
    """
    US at-will employment: termination can happen suddenly with little warning.
    High flexibility for employers; employees face acute but short-lived uncertainty.
    """
    return {
        "employee_stability":  0.35,  # low — sudden layoffs are common
        "company_flexibility": 0.90,  # high — restructuring is fast and legal
        "chronic_stress":      0.65,  # moderate — acute stress, but resolved quickly
        "market_dynamism":     0.85,  # high — active job market, frequent moves
        "career_clarity":      0.70,  # moderate-high — employees learn status quickly
        "economic_efficiency": 0.80,  # high — labor moves to highest-value uses
    }


def score_japan_employment():
    """
    Japan lifetime employment: explicit termination is rare but employees may be
    sidelined, demoted, or pressured — the threat lingers indefinitely.
    """
    return {
        "employee_stability":  0.80,  # high — formal layoffs are uncommon
        "company_flexibility": 0.25,  # low — restructuring is slow and costly
        "chronic_stress":      0.30,  # low score = high chronic stress
        "market_dynamism":     0.30,  # low — mid-career job changes are stigmatized
        "career_clarity":      0.25,  # low — employees rarely receive direct feedback
        "economic_efficiency": 0.45,  # moderate — misallocated workers remain on payroll
    }


def weighted_score(scores):
    """Return the weighted average score across all dimensions using WEIGHTS."""
    return sum(WEIGHTS[dim] * scores[dim] for dim in DIMENSIONS)


def main():
    us = score_us_employment()
    jp = score_japan_employment()

    us_total = weighted_score(us)
    jp_total = weighted_score(jp)

    width = 42
    print("\n" + "=" * width)
    print("      JOB SECURITY CULTURE COMPARISON")
    print("=" * width)
    print(f"{'Dimension':<22} | {'US':>6} | {'Japan':>6}")
    print("-" * width)
    for dim, label in DIMENSIONS.items():
        print(f"{dim:<22} | {us[dim]:6.2f} | {jp[dim]:6.2f}")
    print("-" * width)
    print(f"{'Overall (weighted)':<22} | {us_total:6.2f} | {jp_total:6.2f}")
    print("=" * width)

    delta = us_total - jp_total
    if delta > 0:
        print(
            f"\nConclusion: US-style at-will employment scores higher "
            f"(+{delta:.2f}) on the combined metric.\n"
        )
    elif delta < 0:
        print(
            f"\nConclusion: Japan-style lifetime employment scores higher "
            f"(+{-delta:.2f}) on the combined metric.\n"
        )
    else:
        print("\nConclusion: Both systems score equally overall.\n")

    print(
        "Note: Neither system is universally superior — tradeoffs differ by\n"
        "perspective (employer vs. employee) and socioeconomic context.\n"
    )


if __name__ == "__main__":
    main()
