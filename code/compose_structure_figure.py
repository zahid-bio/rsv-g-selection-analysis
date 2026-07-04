from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
import matplotlib.image as mpimg

PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "figures" / "rsv_g_ccd_structure_pymol_raw.png"
OUT = PROJECT / "figures" / "rsv_g_ccd_structure_pymol.png"

PYMOL_BWR = LinearSegmentedColormap.from_list(
    "pymol_bwr", ["#0000ff", "#ffffff", "#ff0000"])


def main():
    if not RAW.exists():
        raise SystemExit(f"missing {RAW}; run PyMOL render first "
                         f"(pymol -cq code/render_ccd_pymol.pml)")
    img = mpimg.imread(str(RAW))

    fig = plt.figure(figsize=(11, 8.8))
    ax = fig.add_axes([0.02, 0.10, 0.84, 0.82])
    ax.imshow(img)
    ax.axis("off")

    fig.suptitle("RSV-A2 G central conserved domain: the CX3C / cystine-noose core "
                 "is under purifying selection",
                 fontsize=13, fontweight="bold", y=0.975)
    ax.set_title("PDB 5WN9 chain A (residues 169–189), backbone colored by per-codon dN−dS (FEL, RSV-A). "
                 "Gold = disulfides 173–186 & 176–182.",
                 fontsize=9.5, pad=6)

    cax = fig.add_axes([0.88, 0.28, 0.025, 0.46])
    norm = Normalize(vmin=-3, vmax=3)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=PYMOL_BWR), cax=cax)
    cb.set_ticks([-3, -2, -1, 0, 1, 2, 3])
    cb.set_label("per-codon dN − dS  (FEL β − α)", fontsize=10)
    cax.text(0.5, 1.03, "diversifying", transform=cax.transAxes,
             ha="center", va="bottom", fontsize=8.5, color="#b00000")
    cax.text(0.5, -0.03, "purifying", transform=cax.transAxes,
             ha="center", va="top", fontsize=8.5, color="#0000b0")

    fig.text(0.03, 0.045,
             "Peptide spans the CX3C motif and cystine noose only; the full CCR (157–198) and the "
             "intrinsically-disordered mucin domains are not part of any experimental G structure. "
             "dN−dS from this study (rsv_a FEL); colorbar limits match the render (white = 0).",
             fontsize=7.8, color="0.3", wrap=True)

    fig.savefig(OUT, dpi=200)
    print(f"wrote -> {OUT}")


if __name__ == "__main__":
    main()
