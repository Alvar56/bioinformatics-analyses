"""Calculations behind routine food-microbiology results.

* colony_count  - colonies per gram (or ml) from plate counts on successive dilutions,
                  using the weighted mean of ISO 7218.
* mpn           - Most Probable Number from a multiple-tube series, by maximum
                  likelihood, with a likelihood-ratio 95 % confidence interval.

Run as a script for a worked example of each:

    python 04-food-micro-qc/microqc.py
"""
from __future__ import annotations

import math
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Plate counts (ISO 7218, weighted mean over two successive dilutions)
# ---------------------------------------------------------------------------

@dataclass
class Plate:
    dilution: float      # e.g. 1e-2 for the 10^-2 dilution of the sample
    colonies: int


def colony_count(plates: list[Plate], volume_ml: float = 1.0,
                 lower: int = 10, upper: int = 300) -> float:
    """CFU per gram or ml of the original sample.

    N = sum(C) / (V * (n1 + 0.1 * n2) * d)

    where sum(C) is the total of colonies on the retained plates of two successive
    dilutions, n1 and n2 the number of retained plates at the first (less diluted) and
    second dilution, d the first retained dilution and V the inoculated volume.
    Plates are retained when they hold between `lower` and `upper` colonies.
    """
    retained = [p for p in plates if lower <= p.colonies <= upper]
    if not retained:
        raise ValueError(f"no plate between {lower} and {upper} colonies; "
                         "report as an estimate or below the limit")
    levels = sorted({p.dilution for p in retained}, reverse=True)[:2]
    first = levels[0]
    n1 = sum(1 for p in retained if p.dilution == first)
    n2 = sum(1 for p in retained if len(levels) > 1 and p.dilution == levels[1])
    total = sum(p.colonies for p in retained if p.dilution in levels)
    return total / (volume_ml * (n1 + 0.1 * n2) * first)


def round_iso(n: float) -> str:
    """Two significant figures, written as a power of ten (ISO 7218 reporting)."""
    exp = math.floor(math.log10(n))
    mant = round(n / 10 ** exp, 1)
    if mant >= 10:
        mant, exp = mant / 10, exp + 1
    return f"{mant:.1f} × 10^{exp}"


# ---------------------------------------------------------------------------
# Most Probable Number (maximum likelihood)
# ---------------------------------------------------------------------------

def _loglik(lam: float, positives, tubes, amounts) -> float:
    ll = 0.0
    for p, n, v in zip(positives, tubes, amounts):
        x = lam * v
        if p:
            ll += p * math.log(-math.expm1(-x))
        ll -= (n - p) * x
    return ll


def _score(lam: float, positives, tubes, amounts) -> float:
    """d logL / d lambda; strictly decreasing in lambda, so it has one root."""
    s = 0.0
    for p, n, v in zip(positives, tubes, amounts):
        x = lam * v
        if p:
            s += p * v / (-math.expm1(-x))
        s -= n * v
    return s


def mpn(positives, tubes, amounts) -> float:
    """Most Probable Number of organisms per unit of `amounts` (e.g. per gram).

    positives : positive tubes at each level, e.g. (3, 2, 0)
    tubes     : tubes per level, e.g. (3, 3, 3)
    amounts   : sample amount per tube at each level, e.g. (0.1, 0.01, 0.001) g
    Returns 0 if every tube is negative and inf if every tube is positive.
    """
    if sum(positives) == 0:
        return 0.0
    if all(p == n for p, n in zip(positives, tubes)):
        return math.inf
    lo, hi = 1e-9, 1e12
    for _ in range(200):                      # bisection on a log scale
        mid = math.sqrt(lo * hi)
        if _score(mid, positives, tubes, amounts) > 0:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def mpn_interval(positives, tubes, amounts, level: float = 3.841):
    """Likelihood-ratio interval: lambdas with 2*(logL(MPN) - logL) <= chi2(1, 0.95)."""
    best = mpn(positives, tubes, amounts)
    if best in (0.0, math.inf):
        raise ValueError("interval is undefined when all tubes agree")
    top = _loglik(best, positives, tubes, amounts)

    def inside(lam):
        return 2 * (top - _loglik(lam, positives, tubes, amounts)) <= level

    def edge(a, b):                           # a inside, b outside
        for _ in range(200):
            m = math.sqrt(a * b)
            a, b = (m, b) if inside(m) else (a, m)
        return math.sqrt(a * b)

    return edge(best, best / 1e6), edge(best, best * 1e6)


def _demo():
    plates = [Plate(1e-2, 168), Plate(1e-2, 215), Plate(1e-3, 14), Plate(1e-3, 25)]
    n = colony_count(plates)
    print("Plate count (ISO 7218 example): 168, 215 at 10^-2; 14, 25 at 10^-3")
    print(f"  N = {n:,.0f} CFU/g  ->  report {round_iso(n)} CFU/g\n")

    series = ((3, 2, 0), (3, 3, 3), (0.1, 0.01, 0.001))
    m = mpn(*series)
    lo, hi = mpn_interval(*series)
    print("MPN, 3 tubes x 0.1 / 0.01 / 0.001 g, positives 3-2-0")
    print(f"  MPN = {m:.0f} per g, 95 % likelihood interval {lo:.0f}-{hi:.0f}")


if __name__ == "__main__":
    _demo()
