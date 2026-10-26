"""How close is the AlphaFold model of human CCR5 to the experimental structure?

Compares the AlphaFold DB model (UniProt P51681) with the crystal structure of CCR5
bound to maraviroc (PDB 4MBS, chain A). The rubredoxin fusion inserted in 4MBS is
removed, residues are matched by sequence, and the two CA traces are superposed with
iterative outlier rejection. Per-residue deviation is then compared with AlphaFold's
own confidence (pLDDT) and with the transmembrane helices annotated in UniProt.

    python 01-ccr5-model-benchmark/benchmark.py
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common.fetch import alphafold_model, pdb_entry, uniprot_entry  # noqa: E402
from common.plot import house_style, GREEN, GREEN_DARK, AMBER, INK  # noqa: E402
from common.structure import ca_trace, match_by_sequence, rmsd_all, superpose  # noqa: E402

ACC = "P51681"


def main():
    house_style()
    xnum, xseq, xxyz, _ = ca_trace(pdb_entry("4mbs"), "A", keep=lambda n: n < 1000)
    anum, aseq, axyz, plddt = ca_trace(alphafold_model(ACC), "A")
    pairs, ident = match_by_sequence(aseq, xseq)
    P, Q = axyz[pairs[:, 0]], xxyz[pairs[:, 1]]
    d, core, core_rmsd = superpose(P, Q)
    full = rmsd_all(P, Q)

    resnum = anum[pairs[:, 0]]
    conf = plddt[pairs[:, 0]]
    tm = [(f["location"]["start"]["value"], f["location"]["end"]["value"])
          for f in uniprot_entry(ACC)["features"] if f["type"] == "Transmembrane"]
    in_tm = np.zeros(len(resnum), bool)
    for a, b in tm:
        in_tm |= (resnum >= a) & (resnum <= b)
    tm_rmsd = rmsd_all(P[in_tm], Q[in_tm])

    # Figure 1: deviation along the sequence, helices shaded, pLDDT on a second axis.
    fig, ax = plt.subplots(figsize=(10, 3.6))
    for a, b in tm:
        ax.axvspan(a, b, color=GREEN, alpha=0.13, lw=0)
    ax.plot(resnum, d, color=GREEN_DARK, lw=1.6, label="CA deviation after superposition")
    ax.set_xlabel("Residue (UniProt numbering)")
    ax.set_ylabel("Deviation (Å)")
    ax.set_ylim(0, max(8, np.percentile(d, 98) * 1.1))
    ax2 = ax.twinx()
    ax2.plot(resnum, conf, color=AMBER, lw=1.2, label="AlphaFold pLDDT")
    ax2.set_ylabel("pLDDT", color=AMBER)
    ax2.set_ylim(0, 100)
    ax2.spines["right"].set_visible(True)
    ax.set_title("CCR5: AlphaFold model vs crystal structure 4MBS · shaded = transmembrane helices",
                 loc="left", fontsize=10.5, color=INK)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, ncol=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "deviation_along_sequence.png", dpi=200)

    # Figure 2: does pLDDT predict where the model is wrong?
    fig, ax = plt.subplots(figsize=(4.6, 4))
    ax.scatter(conf[~in_tm], d[~in_tm], s=14, color=AMBER, alpha=0.75, label="loops and termini")
    ax.scatter(conf[in_tm], d[in_tm], s=14, color=GREEN_DARK, alpha=0.75, label="transmembrane")
    r = np.corrcoef(conf, d)[0, 1]
    ax.set_xlabel("AlphaFold pLDDT")
    ax.set_ylabel("Deviation from 4MBS (Å)")
    ax.set_title(f"Pearson r = {r:.2f}", loc="left", fontsize=10.5, color=INK)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "plddt_vs_deviation.png", dpi=200)

    rows = [
        ("Residues compared", f"{len(pairs)}"),
        ("Sequence identity over matched residues", f"{ident:.1f} %"),
        ("CA RMSD, all matched residues", f"{full:.2f} Å"),
        ("CA RMSD, transmembrane helices only", f"{tm_rmsd:.2f} Å"),
        ("Core within 4 Å after iterative fit", f"{core.sum()} residues, {core_rmsd:.2f} Å"),
        ("Correlation pLDDT vs deviation", f"r = {r:.2f}"),
    ]
    table = "| Metric | Value |\n| --- | --- |\n" + "\n".join(f"| {a} | {b} |" for a, b in rows) + "\n"
    (HERE / "results.md").write_text(table, encoding="utf-8")
    print(table)


if __name__ == "__main__":
    main()
