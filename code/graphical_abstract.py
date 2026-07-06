"""Generate the graphical abstract: a domain map of RSV-A2 G coloured by the
region-level selection result, summarising the paper's main finding.

Output: figures/graphical_abstract.png
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "figures" / "graphical_abstract.png"

# RSV-A2 G coordinate system
G_LEN = 298
REGIONS = [
    ("N-terminal", 1, 156, "#c9ccd1", "diversifying / variable"),
    ("CCR", 157, 198, "#3b6fb0", "purifying constraint"),
    ("mucin-like C-terminal", 199, 298, "#c0392b", "diversification"),
]
CX3C = (182, 186)

fig, ax = plt.subplots(figsize=(9.2, 4.6))
ax.set_xlim(-6, G_LEN + 6)
ax.set_ylim(0, 10.8)
ax.axis("off")

# Title
ax.text(G_LEN / 2, 9.5,
        "RSV G: the CX3C/cystine-noose core is conserved;\nthe mucin-like flanks diversify",
        ha="center", va="center", fontsize=12.5, fontweight="bold", linespacing=1.2)

# Domain bar
bar_y, bar_h = 5.4, 1.5
for name, start, end, color, _ in REGIONS:
    ax.add_patch(Rectangle((start, bar_y), end - start, bar_h,
                           facecolor=color, edgecolor="black", linewidth=1.1))
    ax.text((start + end) / 2, bar_y + bar_h / 2, name, ha="center", va="center",
            fontsize=9, color="white", fontweight="bold")

# CX3C motif marker within CCR
ax.add_patch(Rectangle((CX3C[0], bar_y), CX3C[1] - CX3C[0], bar_h,
                       facecolor="none", edgecolor="gold", linewidth=2.4))
ax.annotate("CX3C motif + cystine noose (182-186)",
            xy=((CX3C[0] + CX3C[1]) / 2, bar_y + bar_h),
            xytext=(232, bar_y + bar_h + 1.15),
            ha="center", va="bottom", fontsize=8.5, color="#8a6d00",
            arrowprops=dict(arrowstyle="->", color="#8a6d00", lw=1.4))

# Coordinate ticks
for x in (1, 157, 199, 298):
    ax.plot([x, x], [bar_y - 0.25, bar_y], color="black", lw=0.8)
    ax.text(x, bar_y - 0.6, str(x), ha="center", va="top", fontsize=7.5)

# FEL/MEME positive sites: all A2-mappable RSV sites fall in the mucin domain
mucin_sites = [218, 237, 247, 266, 274, 284]
for s in mucin_sites:
    ax.plot(s, bar_y - 0.02, marker="v", color="black", markersize=6, clip_on=False)
ax.text(248, bar_y - 1.5,
        "all A2-mappable FEL/MEME-positive sites (screening-level)",
        ha="center", va="top", fontsize=7.8, style="italic")

# Bottom summary boxes
box_y = 1.2
def box(x, w, title, body, color):
    ax.add_patch(FancyBboxPatch((x, box_y), w, 2.2,
                 boxstyle="round,pad=0.12,rounding_size=0.25",
                 facecolor=color, edgecolor="black", linewidth=1.0, alpha=0.16))
    ax.text(x + w / 2, box_y + 1.62, title, ha="center", va="center",
            fontsize=9, fontweight="bold")
    ax.text(x + w / 2, box_y + 0.72, body, ha="center", va="center", fontsize=8)

box(4, 138, "Conserved core (CCR / CX3C)",
    "No enrichment of diversifying sites.\nCCR more purifying than flanks\n(RSV-A Mann-Whitney p = 0.0066).", "#3b6fb0")
box(156, 138, "Mucin-like flanks + ON1",
    "Where FEL/MEME signals concentrate.\nSite calls do not survive FDR;\nregion-level trend is the primary result.", "#c0392b")

# Footer: pipeline + controls
ax.text(G_LEN / 2, 0.35,
        "Pipeline: MAFFT codon alignment  ->  IQ-TREE ML phylogeny  ->  FEL / MEME (HyPhy)  ->  BH-FDR  |  "
        "GARD: no recombination  |  influenza H3N2 HA positive control (antigenic sites, p = 0.0081)",
        ha="center", va="center", fontsize=7.2, color="#333333")

fig.tight_layout()
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
print(f"wrote -> {OUT}")
