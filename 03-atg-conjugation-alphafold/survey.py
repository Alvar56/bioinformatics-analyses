"""How confidently does AlphaFold model the plant ATG8 conjugation machinery?

ATG8 is attached to membrane lipids by a ubiquitin-like cascade: ATG7 (E1) activates
ATG8 and ATG12, ATG10 (E2) conjugates ATG12 to ATG5, and the ATG12-ATG5-ATG16 complex
guides ATG3 (E2) to transfer ATG8 onto phosphatidylethanolamine.

For each Arabidopsis component this script downloads the AlphaFold DB model, reads the
per-residue confidence (pLDDT, stored in the B-factor column) and annotated domains from
UniProt, and summarises which parts are well-folded and which are likely disordered.
That tells you, before any docking or modelling, which regions are worth trusting.

    python 03-atg-conjugation-alphafold/survey.py
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, BoundaryNorm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common.fetch import alphafold_model, uniprot_entry  # noqa: E402
from common.plot import house_style, INK  # noqa: E402
from common.structure import ca_trace  # noqa: E402

# Arabidopsis thaliana, reviewed UniProt entries
PROTEINS = [
    ("ATG7", "Q94CD5", "E1 · activates ATG8 and ATG12"),
    ("ATG10", "Q8VZ52", "E2 · conjugates ATG12 to ATG5"),
    ("ATG12A", "Q8S924", "ubiquitin-like modifier"),
    ("ATG5", "Q9FFI2", "acceptor of ATG12"),
    ("ATG16", "Q6NNP0", "scaffold · coiled coil + WD40"),
    ("ATG3", "Q0WWQ1", "E2 · transfers ATG8 to PE"),
    ("ATG8A", "Q8LEM4", "ubiquitin-like · lipidated"),
]

# AlphaFold DB confidence bands and their standard colours
BANDS = [0, 50, 70, 90, 100]
BAND_COLOURS = ["#ff7d45", "#ffdb13", "#65cbf3", "#0053d6"]
BAND_NAMES = ["very low (<50)", "low (50–70)", "confident (70–90)", "very high (>90)"]


def domains(entry):
    keep = {"Domain", "Repeat", "Coiled coil", "Motif", "Region"}
    out = []
    for f in entry.get("features", []):
        if f["type"] in keep:
            a, b = f["location"]["start"]["value"], f["location"]["end"]["value"]
            out.append((f["type"], f.get("description") or f["type"], a, b))
    return out


def main():
    house_style()
    rows, profiles = [], []
    for name, acc, role in PROTEINS:
        nums, _, _, plddt = ca_trace(alphafold_model(acc), "A")
        entry = uniprot_entry(acc)
        profiles.append((name, role, nums, plddt, domains(entry)))
        frac = np.histogram(plddt, BANDS)[0] / len(plddt) * 100
        rows.append((name, acc, len(plddt), plddt.mean(), frac[0], frac[3]))

    longest = max(len(p[3]) for p in profiles)
    fig, axes = plt.subplots(len(profiles), 1, figsize=(11, 0.78 * len(profiles) + 1.1), sharex=True)
    cmap = ListedColormap(BAND_COLOURS)
    norm = BoundaryNorm(BANDS, cmap.N)
    for ax, (name, role, nums, plddt, doms) in zip(axes, profiles):
        ax.imshow(plddt[None, :], aspect="auto", cmap=cmap, norm=norm,
                  extent=(nums[0] - 0.5, nums[-1] + 0.5, 0, 1), interpolation="nearest")
        for kind, desc, a, b in doms:
            if kind in ("Domain", "Repeat", "Coiled coil"):
                ax.plot([a, b], [-0.18, -0.18], color=INK, lw=2.2, solid_capstyle="butt", clip_on=False)
        ax.set_xlim(0.5, longest + 0.5)
        ax.set_ylim(-0.3, 1)
        ax.set_yticks([])
        ax.tick_params(axis="x", length=0 if ax is not axes[-1] else 3.5)
        ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.text(-0.01, 0.5, name, transform=ax.transAxes, ha="right", va="center",
                fontsize=10.5, fontweight="bold", color=INK)
        ax.text(1.01, 0.5, role, transform=ax.transAxes, ha="left", va="center", fontsize=8.5, color="#4d5a4f")
    axes[-1].set_xlabel("Residue")
    axes[0].set_title("AlphaFold confidence along the Arabidopsis ATG8 conjugation machinery"
                      " · black bars = annotated domains (UniProt)", loc="left", fontsize=10.5, color=INK)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in BAND_COLOURS]
    fig.legend(handles, BAND_NAMES, loc="lower center", ncol=4, frameon=False, fontsize=9,
               title="pLDDT", title_fontsize=9, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0.04, 0.05, 0.86, 1))
    fig.savefig(HERE / "plddt_strips.png", dpi=200)

    table = ("| Protein | UniProt | Length | Mean pLDDT | Very low (<50) | Very high (>90) |\n"
             "| --- | --- | --- | --- | --- | --- |\n")
    table += "\n".join(f"| {n} | [{a}](https://www.uniprot.org/uniprotkb/{a}) | {l} | {m:.1f} | {lo:.0f} % | {hi:.0f} % |"
                       for n, a, l, m, lo, hi in rows) + "\n"
    (HERE / "results.md").write_text(table, encoding="utf-8")
    print(table)


if __name__ == "__main__":
    main()
