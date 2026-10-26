"""One plotting style for every figure in the repository."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GREEN = "#7fb276"
GREEN_DARK = "#2f6b38"
GREEN_LIGHT = "#a9cf96"
AMBER = "#e8a93a"
INK = "#1d2a20"
GRID = "#e4e9e1"


def house_style():
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "#9aa59a",
        "axes.labelcolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "xtick.color": "#4d5a4f",
        "ytick.color": "#4d5a4f",
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "savefig.bbox": "tight",
    })
