"""Benjamini-Hochberg FDR correction for FEL and MEME site-wise p-values.

Reads the raw HyPhy FEL/MEME JSON outputs, computes Benjamini-Hochberg
false-discovery-rate q-values across all codon sites (per dataset and method),
and writes a tidy table of every screening-level site (unadjusted p <= 0.05)
with its adjusted q-value and whether it survives FDR control at q <= 0.05.

Usage:
    python code/apply_fdr.py [--p-cutoff 0.05] [--q-cutoff 0.05] [--out results/summary_fdr_qvalues.csv]

FEL JSON  MLE columns: [alpha, beta, alpha=beta, LRT, p-value, total branch length]
MEME JSON MLE columns: [alpha, beta-, p-, beta+, p+, LRT, p-value, #branches, ...]
"""
import argparse
import csv
import json
from pathlib import Path

DATASETS = ["flu_h3_ha", "rsv_a", "rsv_b"]
FEL_JSON = {k: f"results/{k}_fel.json" for k in DATASETS}
MEME_JSON = {k: f"results/meme/{k}_meme.json" for k in DATASETS}

# p-value column index in the MLE content rows
FEL_P = 4
MEME_P = 6


def benjamini_hochberg(pvals):
    """Return BH-adjusted q-values, preserving input order."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    running_min = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running_min = min(running_min, pvals[i] * n / rank)
        q[i] = min(running_min, 1.0)
    return q


def load_rows(path):
    return json.load(open(path, encoding="utf-8"))["MLE"]["content"]["0"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--p-cutoff", type=float, default=0.05)
    ap.add_argument("--q-cutoff", type=float, default=0.05)
    ap.add_argument("--out", default="results/summary_fdr_qvalues.csv")
    args = ap.parse_args(argv)

    out_rows = []
    counts = []
    for k in DATASETS:
        # FEL (diversifying: beta > alpha)
        fel = load_rows(FEL_JSON[k])
        fel_p = [float(r[FEL_P]) for r in fel]
        fel_q = benjamini_hochberg(fel_p)
        n_fel_p = n_fel_q = 0
        for i, r in enumerate(fel):
            div = float(r[1]) > float(r[0])
            if float(r[FEL_P]) <= args.p_cutoff and div:
                n_fel_p += 1
                survives = fel_q[i] <= args.q_cutoff
                n_fel_q += int(survives)
                out_rows.append({
                    "dataset": k, "method": "FEL", "aln_codon_pos": i + 1,
                    "p_value": f"{fel_p[i]:.6g}", "q_value": f"{fel_q[i]:.6g}",
                    "survives_fdr": "yes" if survives else "no",
                })
        counts.append((k, "FEL", n_fel_p, n_fel_q))

        # MEME (episodic)
        meme = load_rows(MEME_JSON[k])
        meme_p = [float(r[MEME_P]) for r in meme]
        meme_q = benjamini_hochberg(meme_p)
        n_meme_p = n_meme_q = 0
        for i, r in enumerate(meme):
            episodic = float(r[4]) > 0 and float(r[3]) > float(r[0])
            if float(r[MEME_P]) <= args.p_cutoff and episodic:
                n_meme_p += 1
                survives = meme_q[i] <= args.q_cutoff
                n_meme_q += int(survives)
                out_rows.append({
                    "dataset": k, "method": "MEME", "aln_codon_pos": i + 1,
                    "p_value": f"{meme_p[i]:.6g}", "q_value": f"{meme_q[i]:.6g}",
                    "survives_fdr": "yes" if survives else "no",
                })
        counts.append((k, "MEME", n_meme_p, n_meme_q))

    out_rows.sort(key=lambda d: (d["dataset"], d["method"], d["aln_codon_pos"]))
    out_path = Path(args.out)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "dataset", "method", "aln_codon_pos", "p_value", "q_value", "survives_fdr"])
        w.writeheader()
        w.writerows(out_rows)

    print(f"p-cutoff={args.p_cutoff}  q-cutoff(FDR)={args.q_cutoff}")
    print(f"{'dataset':10s} {'method':5s} {'n(p<=cut)':>9s} {'n(q<=cut)':>9s}")
    for k, m, npv, nqv in counts:
        print(f"{k:10s} {m:5s} {npv:9d} {nqv:9d}")
    print(f"wrote -> {out_path}")


if __name__ == "__main__":
    main()
