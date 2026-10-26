"""Figure: the likelihood behind an MPN result, for several tube outcomes.

    python 04-food-micro-qc/plot_mpn.py
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from common.plot import house_style, GREEN, GREEN_DARK, AMBER, INK  # noqa: E402
from microqc import _loglik, mpn, mpn_interval  # noqa: E402

TUBES = (3, 3, 3)
AMOUNTS = (0.1, 0.01, 0.001)
OUTCOMES = [((3, 0, 0), GREEN), ((3, 2, 0), GREEN_DARK), ((3, 3, 1), AMBER)]


def main():
    house_style()
    lam = np.logspace(0, 4, 600)
    fig, ax = plt.subplots(figsize=(8, 3.8))
    for pos, colour in OUTCOMES:
        m = mpn(pos, TUBES, AMOUNTS)
        lo, hi = mpn_interval(pos, TUBES, AMOUNTS)
        top = _loglik(m, pos, TUBES, AMOUNTS)
        rel = np.array([np.exp(_loglik(x, pos, TUBES, AMOUNTS) - top) for x in lam])
        label = f"{'-'.join(map(str, pos))}: MPN {m:.0f}/g  (95 % {lo:.0f}–{hi:.0f})"
        ax.plot(lam, rel, color=colour, lw=2, label=label)
        ax.axvspan(lo, hi, color=colour, alpha=0.08, lw=0)
        ax.axvline(m, color=colour, lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("Organisms per gram")
    ax.set_ylabel("Relative likelihood")
    ax.set_ylim(0, 1.05)
    ax.set_title("Most Probable Number by maximum likelihood · 3 tubes × 0.1 / 0.01 / 0.001 g",
                 loc="left", fontsize=10.5, color=INK)
    ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "mpn_likelihood.png", dpi=200)


if __name__ == "__main__":
    main()
