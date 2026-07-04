import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
RES_DIR = PROJECT / "results"
FIG_DIR = PROJECT / "figures"
FIG_DIR.mkdir(exist_ok=True)

ORDER = ["Nterm", "CCR", "Mucin"]
COLORS = {"Nterm": "#9ecae1", "CCR": "#e6550d", "Mucin": "#a1d99b"}


def load_percodon(label):
    path = RES_DIR / f"{label}_region_selection_percodon.csv"
    by_region = {r: [] for r in ORDER}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            reg = row["region"]
            if reg in by_region:
                by_region[reg].append(float(row["dNdS"]))
    return by_region


def panel(ax, label, title):
    by_region = load_percodon(label)
    tests = json.load(open(RES_DIR / f"{label}_region_selection.json"))["tests"]

    data = [by_region[r] for r in ORDER]
    positions = range(len(ORDER))
    bp = ax.boxplot(data, positions=list(positions), widths=0.55,
                    showfliers=False, patch_artist=True, medianprops=dict(color="black"))
    for patch, reg in zip(bp["boxes"], ORDER):
        patch.set_facecolor(COLORS[reg])
        patch.set_alpha(0.55 if reg == "CCR" else 0.35)
    for i, reg in enumerate(ORDER):
        y = np.asarray(data[i], float)
        x = np.full_like(y, i) + (np.linspace(-0.16, 0.16, len(y)) if len(y) else 0)
        ax.scatter(x, y, s=7, color=COLORS[reg], edgecolor="none", alpha=0.6, zorder=3)

    ax.axhline(0, color="0.5", lw=0.8, ls="--", zorder=1)
    ax.set_xticks(list(positions))
    ax.set_xticklabels([f"{r}\n(n={len(by_region[r])})" for r in ORDER])
    ax.set_ylabel("dN − dS  (FEL β − α)")
    ax.set_title(title, fontweight="bold")

    kw = tests.get("kruskal_Nterm_CCR_Mucin", {})
    cm = tests.get("CCR_vs_Mucin_dNdS_less", {})
    lines = []
    if kw:
        lines.append(f"Kruskal-Wallis: H={kw['H']:.2f}, p={kw['p']:.3g}")
    if cm:
        lines.append(f"CCR<Mucin (MWU): p={cm['p_one_sided_CCR_lower']:.3g}, "
                     f"rbc={cm['rank_biserial']:+.2f}")
    ax.text(0.02, 0.02, "\n".join(lines), transform=ax.transAxes,
            fontsize=8, va="bottom", ha="left",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    panel(axes[0], "rsv_a", "RSV-A  (per-codon selection by region)")
    panel(axes[1], "rsv_b", "RSV-B  (per-codon selection by region)")
    fig.suptitle("Region-level selection pressure on RSV G: CCR is more constrained than the mucin domain",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG_DIR / "rsv_region_selection_dnds.png"
    fig.savefig(out, dpi=200)
    print(f"wrote -> {out}")


if __name__ == "__main__":
    main()
