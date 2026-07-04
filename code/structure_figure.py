import csv
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from Bio.PDB import PDBParser

PROJECT = Path(__file__).resolve().parents[1]
PDB = PROJECT / "data" / "5WN9.pdb"
CSV = PROJECT / "results" / "rsv_a_region_selection_percodon.csv"
OUT = PROJECT / "figures" / "rsv_g_ccd_structure_dnds.png"
CX3C = range(182, 187)
DISULFIDES = [(173, 186), (176, 182)]


def load_dnds():
    d = {}
    for row in csv.DictReader(open(CSV)):
        if row["a2"]:
            d[int(row["a2"])] = float(row["dNdS"])
    return d


def main():
    dnds = load_dnds()
    chain = PDBParser(QUIET=True).get_structure("g", str(PDB))[0]["A"]
    res = [r for r in chain if r.id[0] == " " and "CA" in r]
    nums = [r.id[1] for r in res]
    ca = np.array([r["CA"].coord for r in res])
    vals = np.array([dnds.get(n, 0.0) for n in nums])
    sg = {r.id[1]: r["SG"].coord for r in chain if r.resname == "CYS" and "SG" in r}

    num_to_ca = {n: ca[i] for i, n in enumerate(nums)}
    norm = TwoSlopeNorm(vmin=-5, vcenter=0.0, vmax=1)
    cmap = plt.get_cmap("RdBu_r")

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot(ca[:, 0], ca[:, 1], ca[:, 2], color="0.6", lw=1.5, zorder=1)
    sizes = [130 if n in CX3C else 70 for n in nums]
    edges = ["black" if n in CX3C else "0.4" for n in nums]
    ax.scatter(ca[:, 0], ca[:, 1], ca[:, 2], c=vals, cmap=cmap, norm=norm,
               s=sizes, edgecolors=edges, linewidths=[1.6 if n in CX3C else 0.6 for n in nums],
               depthshade=False, zorder=3)

    for a, b in DISULFIDES:
        if a in num_to_ca and b in num_to_ca:
            p, q = num_to_ca[a], num_to_ca[b]
            ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]],
                    color="#e6b800", lw=3.0, ls=(0, (1, 1)), zorder=2)

    ax.text(*ca[0], f"  N ({nums[0]})", fontsize=9, color="black")
    ax.text(*ca[-1], f"  C ({nums[-1]})", fontsize=9, color="black")
    cx = np.array([ca[i] for i, n in enumerate(nums) if n in CX3C]).mean(axis=0)
    ax.text(*cx, "  CX3C\n  (182-186)", fontsize=10, fontweight="bold", color="black")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("per-codon dN − dS (FEL)\n← purifying      diversifying →")

    ax.set_title("RSV-A2 G central conserved domain (PDB 5WN9, res 169–189)\n"
                 "backbone colored by selection; yellow = cystine-noose disulfides",
                 fontsize=11, fontweight="bold")
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    ax.view_init(elev=18, azim=-60)
    ax.text2D(0.02, 0.02,
              "Large black-outlined spheres = CX3C motif. Numbering = RSV-A2 (M11486).\n"
              "Structure colored by dN−dS from this study's FEL analysis (rsv_a).",
              transform=ax.transAxes, fontsize=8, va="bottom",
              bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.9))

    fig.tight_layout()
    fig.savefig(OUT, dpi=200)
    print(f"wrote -> {OUT}")

    from Bio.PDB import PDBIO, Select
    for r in chain:
        v = dnds.get(r.id[1], 0.0)
        for atom in r:
            atom.set_bfactor(round(v, 2))

    class ChainA(Select):
        def accept_chain(self, ch):
            return ch.id == "A"

    io = PDBIO()
    io.set_structure(chain.get_parent().get_parent())
    bpath = PROJECT / "data" / "5WN9_chainA_dnds_bfactor.pdb"
    io.save(str(bpath), ChainA())
    print(f"wrote -> {bpath}")


if __name__ == "__main__":
    main()
