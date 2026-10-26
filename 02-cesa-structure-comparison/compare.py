"""Can structures from other plants stand in for Arabidopsis CESA6?

There is no experimental structure of Arabidopsis CESA6. Structure-based work on it has
to borrow one: an experimental structure of a related cellulose synthase, or a
prediction. This script measures how interchangeable those references are, pairwise:

  * AlphaFold model of Arabidopsis CESA6 (UniProt Q94JQ6)
  * crystal structure of the Arabidopsis CESA3 catalytic domain (PDB 7CK3, chain A)
  * cryo-EM structure of soybean CesA1 (PDB 8VHZ, chain A)
  * cryo-EM structure of poplar CesA8 (PDB 6WLB, chain A)

For each pair, residues are matched by global sequence alignment and the CA traces are
superposed with iterative outlier rejection (a simple, transparent stand-in for
TM-align or FATCAT). The catalytic glycosyltransferase core is compared separately
using the residues around the conserved QxxRW motif.

    python 02-cesa-structure-comparison/compare.py
"""
import itertools
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common.fetch import alphafold_model, pdb_entry  # noqa: E402
from common.plot import house_style, GREEN_DARK, INK  # noqa: E402
from common.structure import ca_trace, match_by_sequence, superpose  # noqa: E402

MODELS = {
    "CESA6 · AlphaFold\n(Arabidopsis)": lambda: ca_trace(alphafold_model("Q94JQ6"), "A"),
    "CESA3 · 7CK3\n(Arabidopsis, X-ray)": lambda: ca_trace(pdb_entry("7ck3"), "A"),
    "CesA1 · 8VHZ\n(soybean, cryo-EM)": lambda: ca_trace(pdb_entry("8vhz"), "A"),
    "CesA8 · 6WLB\n(poplar, cryo-EM)": lambda: ca_trace(pdb_entry("6wlb"), "A"),
}


def main():
    house_style()
    names = list(MODELS)
    data = {n: MODELS[n]() for n in names}
    k = len(names)
    rmsd = np.full((k, k), np.nan)
    ident = np.full((k, k), np.nan)
    rows = []
    for i, j in itertools.combinations(range(k), 2):
        _, sa, xa, _ = data[names[i]]
        _, sb, xb, _ = data[names[j]]
        pairs, pid = match_by_sequence(sa, sb)
        d, core, r = superpose(xa[pairs[:, 0]], xb[pairs[:, 1]])
        rmsd[i, j] = rmsd[j, i] = r
        ident[i, j] = ident[j, i] = pid
        cov = core.sum() / min(len(sa), len(sb)) * 100
        rows.append((names[i].replace("\n", " "), names[j].replace("\n", " "),
                     len(pairs), pid, core.sum(), r, cov))

    # Where is the QxxRW motif (catalytic core landmark) in each chain?
    motif = {n: [m.start() for m in re.finditer(r"Q..RW", data[n][1])] for n in names}

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    labels = [n.split("\n")[0] for n in names]
    for ax, mat, title, fmt, cmap, lo, hi in (
        (axes[0], rmsd, "Core CA RMSD after superposition (Å)", "{:.2f}", "Greens", 0.8, 1.4),
        (axes[1], ident, "Sequence identity over aligned residues (%)", "{:.0f}", "YlOrBr", 60, 95),
    ):
        shown = np.ma.masked_invalid(mat)
        im = ax.imshow(shown, cmap=cmap, vmin=lo, vmax=hi)
        ax.set_xticks(range(k), labels, rotation=30, ha="right")
        ax.set_yticks(range(k), labels)
        ax.grid(False)
        for a in range(k):
            for b in range(k):
                txt = "—" if a == b else fmt.format(mat[a, b])
                ax.text(b, a, txt, ha="center", va="center", fontsize=10,
                        color="white" if a != b and (mat[a, b] - lo) / (hi - lo) > 0.6 else INK)
        ax.set_title(title, loc="left", fontsize=10.5, color=INK)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(HERE / "pairwise_matrix.png", dpi=200)

    table = ("| Structure A | Structure B | Aligned residues | Identity | Core residues | Core RMSD | Core coverage |\n"
             "| --- | --- | --- | --- | --- | --- | --- |\n")
    table += "\n".join(f"| {a} | {b} | {n} | {p:.0f} % | {c} | {r:.2f} Å | {cv:.0f} % |"
                       for a, b, n, p, c, r, cv in rows) + "\n"
    table += "\nQxxRW motif found at (0-based index in the modelled chain): " + \
        "; ".join(f"{n.split(chr(10))[0]}: {v}" for n, v in motif.items()) + "\n"
    (HERE / "results.md").write_text(table, encoding="utf-8")
    print(table)


if __name__ == "__main__":
    main()
