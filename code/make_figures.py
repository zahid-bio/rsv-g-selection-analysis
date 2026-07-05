import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
import functional_regions as fr

PROJECT = Path(__file__).resolve().parents[1]
FIG_DIR = PROJECT / "figures"
FIG_DIR.mkdir(exist_ok=True)
RES_DIR = PROJECT / "results"


def null_histogram(enrichment_json_path, out_png, title, region_label):
    d = json.load(open(enrichment_json_path))
    null = d["null_distribution"]
    obs = d["observed_fraction"]
    n = d["n_selected_in_domain"]
    exp = d["expected_fraction"]
    p = d["p_value"]

    fig, ax = plt.subplots(figsize=(6.0, 4.0), dpi=150)
    if null:
        ax.hist([x * 100 for x in null], bins=30, color="#8ea9d8", edgecolor="#3b64a6",
                alpha=0.85, label="null (random placement)")
    ax.axvline(obs * 100, color="crimson", linewidth=2.2,
               label=f"observed = {obs*100:.1f}%  (n={n})")
    ax.axvline(exp * 100, color="black", linestyle="--", linewidth=1.2,
               label=f"expected = {exp*100:.1f}%")
    ax.set_xlabel(f"Fraction of selected sites in {region_label}  (%)")
    ax.set_ylabel("Number of permutations")
    ax.set_title(f"{title}\np = {p:.4f}   (10,000 permutations, seed=42)")
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    print(f"  wrote {out_png}")


def rsv_lollipop(out_png):
    def read_sites(path, source):
        rows = []
        for r in csv.DictReader(open(path)):
            a2 = r.get("a2_residue", "")
            if a2 == "":
                continue
            rows.append((int(a2), float(r["beta_dN"]), float(r["p_value"]), source))
        return rows

    A = read_sites(RES_DIR / "rsv_a_selected_a2.csv", "A")
    B = read_sites(RES_DIR / "rsv_b_selected_a2.csv", "B")

    fig, ax = plt.subplots(figsize=(9.0, 3.6), dpi=150)
    ax.axhline(0, color="#666", linewidth=1.0)
    ax.axvspan(min(fr.RSV_G_CCR), max(fr.RSV_G_CCR),
               ymin=0.15, ymax=0.85, color="#c6dbef", alpha=0.55,
               label=f"CCR (A2 {min(fr.RSV_G_CCR)}-{max(fr.RSV_G_CCR)})")
    ax.axvspan(min(fr.RSV_G_CX3C), max(fr.RSV_G_CX3C),
               ymin=0.15, ymax=0.85, color="#fdae6b", alpha=0.85,
               label=f"CX3C (A2 {min(fr.RSV_G_CX3C)}-{max(fr.RSV_G_CX3C)})")
    for a2, dN, p, s in A:
        ax.vlines(a2, 0, dN, color="#3182bd", linewidth=1.5)
        ax.plot([a2], [dN], "o", color="#3182bd", markersize=6, zorder=3)
        ax.annotate(str(a2), (a2, dN), textcoords="offset points",
                    xytext=(0, 5), ha="center", fontsize=7, color="#08306b")
    for a2, dN, p, s in B:
        ax.vlines(a2, 0, -dN, color="#e6550d", linewidth=1.5)
        ax.plot([a2], [-dN], "s", color="#e6550d", markersize=6, zorder=3)
        ax.annotate(str(a2), (a2, -dN), textcoords="offset points",
                    xytext=(0, -12), ha="center", fontsize=7, color="#8c2d04")

    ax.set_xlim(0, fr.RSV_G_LENGTH + 5)
    max_dN = max([d for _, d, *_ in A + B] + [1.0])
    ax.set_ylim(-max_dN - 1, max_dN + 1)
    ax.set_xlabel("RSV G residue position (A2 numbering)")
    ax.set_ylabel(r"$\beta$ (dN)  — RSV-A up, RSV-B down")
    ax.set_title("Positively-selected sites in RSV G (FEL p <= 0.05)\n"
                 "compared to pre-specified CCR and CX3C motif")
    a_patch = mpatches.Patch(color="#3182bd", label="RSV-A")
    b_patch = mpatches.Patch(color="#e6550d", label="RSV-B")
    ccr_patch = mpatches.Patch(color="#c6dbef", label="CCR (pre-spec.)", alpha=0.55)
    cx3c_patch = mpatches.Patch(color="#fdae6b", label="CX3C (pre-spec.)", alpha=0.85)
    ax.legend(handles=[a_patch, b_patch, ccr_patch, cx3c_patch], loc="upper left",
              fontsize=8, frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    print(f"  wrote {out_png}")


def tree_figure(treefile, out_png, title):
    tree = Phylo.read(treefile, "newick")
    try:
        tree.ladderize()
    except Exception:
        pass
    n_taxa = tree.count_terminals()
    fig_h = max(4.0, min(0.14 * n_taxa, 18.0))
    fig, ax = plt.subplots(figsize=(8.0, fig_h), dpi=150)
    Phylo.draw(tree, do_show=False, axes=ax, show_confidence=False,
               label_func=lambda c: c.name if c.is_terminal() else "")
    ax.set_title(f"{title}  (n={n_taxa} taxa)")
    ax.tick_params(labelsize=6)
    for lbl in ax.get_yticklabels():
        lbl.set_fontsize(5)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    print(f"  wrote {out_png}  (taxa={n_taxa})")


def summary_tables():
    per_site = []
    for r in csv.DictReader(open(RES_DIR / "flu_h3_ha_selected_h3.csv")):
        h3 = r.get("h3_ha1_number", "")
        if h3 in ("", "None"):
            continue
        h3 = int(h3)
        site = ""
        for k, res in fr.H3_ANTIGENIC_SITES.items():
            if h3 in res:
                site = k
                break
        per_site.append({
            "target": "flu (H3N2)", "aln_pos": r["aln_codon_pos"],
            "reference_residue": h3, "reference": "H3 mature HA1",
            "region_hit": f"Site {site}" if site else "outside",
            "beta_dN": r["beta_dN"], "p_value": r["p_value"],
        })
    for r in csv.DictReader(open(RES_DIR / "rsv_a_selected_a2.csv")):
        a2 = r.get("a2_residue", "")
        if a2 == "":
            per_site.append({
                "target": "RSV-A", "aln_pos": r["aln_codon_pos"],
                "reference_residue": "n/a (ON1 insertion)",
                "reference": "RSV-A2 (M11486)",
                "region_hit": "n/a", "beta_dN": r["beta_dN"], "p_value": r["p_value"],
            })
            continue
        a2 = int(a2)
        rhit = "CX3C" if a2 in fr.RSV_G_CX3C else ("CCR" if a2 in fr.RSV_G_CCR else "outside")
        per_site.append({
            "target": "RSV-A", "aln_pos": r["aln_codon_pos"],
            "reference_residue": a2, "reference": "RSV-A2 (M11486)",
            "region_hit": rhit, "beta_dN": r["beta_dN"], "p_value": r["p_value"],
        })
    for r in csv.DictReader(open(RES_DIR / "rsv_b_selected_a2.csv")):
        a2 = r.get("a2_residue", "")
        if a2 == "":
            per_site.append({
                "target": "RSV-B", "aln_pos": r["aln_codon_pos"],
                "reference_residue": "n/a", "reference": "RSV-A2 (M11486)",
                "region_hit": "n/a", "beta_dN": r["beta_dN"], "p_value": r["p_value"],
            })
            continue
        a2 = int(a2)
        rhit = "CX3C" if a2 in fr.RSV_G_CX3C else ("CCR" if a2 in fr.RSV_G_CCR else "outside")
        per_site.append({
            "target": "RSV-B", "aln_pos": r["aln_codon_pos"],
            "reference_residue": a2, "reference": "RSV-A2 (M11486)",
            "region_hit": rhit, "beta_dN": r["beta_dN"], "p_value": r["p_value"],
        })

    site_csv = RES_DIR / "summary_selected_sites.csv"
    with open(site_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["target", "aln_pos", "reference_residue",
                                          "reference", "region_hit", "beta_dN", "p_value"])
        w.writeheader()
        for row in per_site:
            w.writerow(row)
    print(f"  wrote {site_csv}")

    per_target = []
    for target, path, region_label in [
        ("flu",             "flu_enrichment.json",             "H3 antigenic sites A-E"),
        ("rsv_a (CCR)",     "rsv_a_enrichment.json",           "RSV G CCR 157-198"),
        ("rsv_a (CX3C)",    "rsv_a_enrichment_cx3c.json",      "RSV G CX3C 182-186"),
        ("rsv_b (CCR)",     "rsv_b_enrichment.json",           "RSV G CCR 157-198"),
        ("rsv_b (CX3C)",    "rsv_b_enrichment_cx3c.json",      "RSV G CX3C 182-186"),
        ("rsv_combined (CCR)",  "rsv_combined_enrichment.json",     "RSV G CCR 157-198"),
        ("rsv_combined (CX3C)", "rsv_combined_enrichment_cx3c.json","RSV G CX3C 182-186"),
    ]:
        j = json.load(open(RES_DIR / path))
        per_target.append({
            "target": target,
            "region": region_label,
            "n_sites_tested": j["n_selected_in_domain"],
            "hits": j["hits_in_region"],
            "observed_frac_pct": round(j["observed_fraction"] * 100, 2),
            "expected_frac_pct": round(j["expected_fraction"] * 100, 2),
            "enrichment_ratio": round(j["enrichment_ratio"], 2)
                if j["enrichment_ratio"] != float("inf") else "inf",
            "p_value": round(j["p_value"], 4),
            "significant_at_0.05": "YES" if j["p_value"] < 0.05 else "no",
        })
    enr_csv = RES_DIR / "summary_enrichment.csv"
    with open(enr_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(per_target[0].keys()))
        w.writeheader()
        for row in per_target:
            w.writerow(row)
    print(f"  wrote {enr_csv}")

    md = RES_DIR / "summary_enrichment.md"
    with open(md, "w", encoding="utf-8") as f:
        f.write("# Enrichment summary — U6 project\n\n")
        f.write("| Target | Region | n | hits | obs % | exp % | enrichment | p-value | sig |\n")
        f.write("|---|---|---:|---:|---:|---:|---:|---:|:---:|\n")
        for r in per_target:
            f.write(f"| {r['target']} | {r['region']} | {r['n_sites_tested']} | "
                    f"{r['hits']} | {r['observed_frac_pct']} | {r['expected_frac_pct']} | "
                    f"{r['enrichment_ratio']}x | {r['p_value']} | {r['significant_at_0.05']} |\n")
        f.write("\nAll tests: 10,000 permutations, seed 42, one-sided (enrichment).\n")

    region_rows = []
    for target, path in [("rsv_a", "rsv_a_region_selection.json"),
                         ("rsv_b", "rsv_b_region_selection.json")]:
        rp = RES_DIR / path
        if not rp.exists():
            continue
        t = json.load(open(rp))["tests"]
        kw = t.get("kruskal_Nterm_CCR_Mucin", {})
        cm = t.get("CCR_vs_Mucin_dNdS_less", {})
        cn = t.get("CCR_vs_nonCCR_dNdS_less", {})
        region_rows.append({
            "target": target,
            "comparison": "Nterm/CCR/Mucin (omnibus)",
            "test": "Kruskal-Wallis",
            "statistic": f"H={kw.get('H'):.2f}" if kw else "",
            "effect_size": f"eps2={kw.get('epsilon_squared'):.3f}" if kw else "",
            "p_value": round(kw["p"], 4) if kw else "",
            "significant_at_0.05": ("YES" if kw and kw["p"] < 0.05 else "no"),
        })
        region_rows.append({
            "target": target,
            "comparison": "CCR < Mucin (more purifying)",
            "test": "Mann-Whitney U (1-sided)",
            "statistic": f"U={cm.get('U'):.0f}" if cm and cm.get("U") is not None else "",
            "effect_size": f"rbc={cm.get('rank_biserial'):+.3f}" if cm else "",
            "p_value": round(cm["p_one_sided_CCR_lower"], 4) if cm else "",
            "significant_at_0.05": ("YES" if cm and cm["p_one_sided_CCR_lower"] < 0.05 else "no"),
        })
        region_rows.append({
            "target": target,
            "comparison": "CCR < non-CCR (mapped)",
            "test": "Mann-Whitney U (1-sided)",
            "statistic": f"U={cn.get('U'):.0f}" if cn and cn.get("U") is not None else "",
            "effect_size": f"rbc={cn.get('rank_biserial'):+.3f}" if cn else "",
            "p_value": round(cn["p_one_sided_CCR_lower"], 4) if cn else "",
            "significant_at_0.05": ("YES" if cn and cn["p_one_sided_CCR_lower"] < 0.05 else "no"),
        })
    if region_rows:
        reg_csv = RES_DIR / "summary_region_selection.csv"
        with open(reg_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(region_rows[0].keys()))
            w.writeheader()
            for row in region_rows:
                w.writerow(row)
        print(f"  wrote {reg_csv}")
        with open(md, "a", encoding="utf-8") as f:
            f.write("\n## Region-level selection (threshold-free companion)\n\n")
            f.write("Per-codon dN-dS (FEL beta-alpha) compared across regions using ALL codons "
                    "(not just the FEL-significant handful). Negative effect size = CCR "
                    "under stronger purifying selection.\n\n")
            f.write("| Target | Comparison | Test | Statistic | Effect size | p-value | sig |\n")
            f.write("|---|---|---|---|---|---:|:---:|\n")
            for r in region_rows:
                f.write(f"| {r['target']} | {r['comparison']} | {r['test']} | "
                        f"{r['statistic']} | {r['effect_size']} | {r['p_value']} | "
                        f"{r['significant_at_0.05']} |\n")
            f.write("\nRank-biserial (rbc): 0.1 small, 0.3 medium, 0.5 large. "
                    "RSV-B A2 mapping loses 4 C-terminal mucin residues "
                    "(a2_ref_len 294 vs 298; MAFFT --keeplength).\n")
    print(f"  wrote {md}")


def main():
    print("[1/4] permutation-null histograms")
    null_histogram(RES_DIR / "flu_enrichment.json",
                   FIG_DIR / "flu_null_hist.png",
                   "Flu H3 HA — selected sites vs. random placement",
                   "H3 antigenic sites A-E")
    null_histogram(RES_DIR / "rsv_a_enrichment.json",
                   FIG_DIR / "rsv_a_ccr_null_hist.png",
                   "RSV-A G — selected sites vs. random placement",
                   "CCR (A2 157-198)")
    null_histogram(RES_DIR / "rsv_a_enrichment_cx3c.json",
                   FIG_DIR / "rsv_a_cx3c_null_hist.png",
                   "RSV-A G — selected sites vs. random placement",
                   "CX3C motif (A2 182-186)")
    null_histogram(RES_DIR / "rsv_b_enrichment.json",
                   FIG_DIR / "rsv_b_ccr_null_hist.png",
                   "RSV-B G — selected sites vs. random placement",
                   "CCR (A2 157-198)")
    null_histogram(RES_DIR / "rsv_b_enrichment_cx3c.json",
                   FIG_DIR / "rsv_b_cx3c_null_hist.png",
                   "RSV-B G — selected sites vs. random placement",
                   "CX3C motif (A2 182-186)")
    null_histogram(RES_DIR / "rsv_combined_enrichment.json",
                   FIG_DIR / "rsv_combined_ccr_null_hist.png",
                   "RSV-A+B G (pooled) — selected sites vs. random placement",
                   "CCR (A2 157-198)")
    null_histogram(RES_DIR / "rsv_combined_enrichment_cx3c.json",
                   FIG_DIR / "rsv_combined_cx3c_null_hist.png",
                   "RSV-A+B G (pooled) — selected sites vs. random placement",
                   "CX3C motif (A2 182-186)")

    print("[2/4] RSV G lollipop plot")
    rsv_lollipop(FIG_DIR / "rsv_selected_sites_lollipop.png")

    print("[3/4] phylogenetic trees")
    tree_figure(RES_DIR / "flu_h3_ha_codon_iqtree.treefile",
                FIG_DIR / "flu_h3_ha_tree.png",
                "Flu H3N2 HA phylogeny (IQ-TREE, TVM+F+I+R2)")
    tree_figure(RES_DIR / "rsv_a_codon_iqtree.treefile",
                FIG_DIR / "rsv_a_tree.png",
                "RSV-A G phylogeny (IQ-TREE, TN+F+I+R2)")
    tree_figure(RES_DIR / "rsv_b_codon_iqtree.treefile",
                FIG_DIR / "rsv_b_tree.png",
                "RSV-B G phylogeny (IQ-TREE, TN+F+I+R2)")

    print("[4/4] summary tables")
    summary_tables()

    print("\nAll figures + tables written under figures/ and results/.")


if __name__ == "__main__":
    main()
