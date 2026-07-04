import sys, json, argparse
from pathlib import Path
from Bio import SeqIO
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import functional_regions as fr

CCR = set(fr.RSV_G_CCR)
CX3C = set(fr.RSV_G_CX3C)
NTERM = set(range(1, 157))
MUCIN = set(range(199, fr.RSV_G_LENGTH + 1))

PCUT = 0.05


def build_aln_to_a2(aln_path, ref_id):
    ref = next((r for r in SeqIO.parse(aln_path, "fasta") if r.id == ref_id), None)
    if ref is None:
        raise SystemExit(f"reference id {ref_id!r} not found in {aln_path}")
    mapping, a2 = {}, 0
    for col, aa in enumerate(str(ref.seq), start=1):
        if aa == "-":
            mapping[col] = None
        else:
            a2 += 1
            mapping[col] = a2
    return mapping, a2


def region_of(a2):
    if a2 is None:
        return "insertion"
    if a2 in CCR:
        return "CCR"
    if a2 in MUCIN:
        return "Mucin"
    if a2 in NTERM:
        return "Nterm"
    return "other"


def summarise(dnds, betas, alphas, ps):
    dnds = np.asarray(dnds, float)
    diversifying = int(np.sum((np.asarray(betas) > np.asarray(alphas)) & (np.asarray(ps) <= PCUT)))
    purifying = int(np.sum((np.asarray(betas) < np.asarray(alphas)) & (np.asarray(ps) <= PCUT)))
    n = len(dnds)
    return {
        "n_codons": n,
        "mean_dN_minus_dS": float(np.mean(dnds)) if n else None,
        "median_dN_minus_dS": float(np.median(dnds)) if n else None,
        "n_diversifying_p05": diversifying,
        "n_purifying_p05": purifying,
        "frac_purifying_p05": (purifying / n) if n else None,
        "frac_diversifying_p05": (diversifying / n) if n else None,
    }


def load_fel_by_region(fel_json, aln_to_a2):
    d = json.load(open(fel_json))
    rows = d["MLE"]["content"]["0"]
    buckets = {}
    per_codon = []
    for i, r in enumerate(rows, start=1):
        alpha, beta = r[0], r[1]
        p = r[4]
        a2 = aln_to_a2.get(i)
        reg = region_of(a2)
        rec = {"aln_col": i, "a2": a2, "region": reg,
               "alpha": alpha, "beta": beta, "dNdS": beta - alpha, "p": p}
        per_codon.append(rec)
        b = buckets.setdefault(reg, {"dnds": [], "beta": [], "alpha": [], "p": []})
        b["dnds"].append(beta - alpha)
        b["beta"].append(beta)
        b["alpha"].append(alpha)
        b["p"].append(p)
    return per_codon, buckets


def mwu(a, b, alternative):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return None, None, None
    U, p = stats.mannwhitneyu(a, b, alternative=alternative)
    rbc = 2.0 * U / (len(a) * len(b)) - 1.0
    return float(U), float(p), float(rbc)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fel-json", required=True)
    ap.add_argument("--aln-with-A2", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--ref-id", default="M11486_A2")
    args = ap.parse_args(argv)

    aln_to_a2, ref_len = build_aln_to_a2(args.aln_with_A2, args.ref_id)
    per_codon, buckets = load_fel_by_region(args.fel_json, aln_to_a2)

    regions = {}
    for reg, b in buckets.items():
        regions[reg] = summarise(b["dnds"], b["beta"], b["alpha"], b["p"])

    cx3c_codons = [c for c in per_codon if c["a2"] in CX3C]
    if cx3c_codons:
        regions["CX3C_nested"] = summarise(
            [c["dNdS"] for c in cx3c_codons],
            [c["beta"] for c in cx3c_codons],
            [c["alpha"] for c in cx3c_codons],
            [c["p"] for c in cx3c_codons],
        )

    dnds = {reg: [c["dNdS"] for c in per_codon if c["region"] == reg] for reg in
            ("Nterm", "CCR", "Mucin", "insertion")}
    ccr = dnds["CCR"]
    mucin = dnds["Mucin"]
    nonccr = [c["dNdS"] for c in per_codon if c["region"] != "CCR" and c["a2"] is not None]
    cx3c_dnds = [c["dNdS"] for c in cx3c_codons]
    ccr_no_cx3c = [c["dNdS"] for c in per_codon if c["region"] == "CCR" and c["a2"] not in CX3C]

    tests = {}
    kw_groups = [dnds["Nterm"], ccr, mucin]
    kw_labels = ["Nterm", "CCR", "Mucin"]
    if all(len(g) > 0 for g in kw_groups):
        H, p_kw = stats.kruskal(*kw_groups)
        n_tot = sum(len(g) for g in kw_groups)
        eps2 = float((H - len(kw_groups) + 1) / (n_tot - len(kw_groups))) if n_tot > len(kw_groups) else None
        tests["kruskal_Nterm_CCR_Mucin"] = {
            "H": float(H), "p": float(p_kw), "df": len(kw_groups) - 1,
            "epsilon_squared": eps2, "groups": kw_labels,
            "n_per_group": {lab: len(g) for lab, g in zip(kw_labels, kw_groups)},
        }

    U, p, rbc = mwu(ccr, mucin, "less")
    tests["CCR_vs_Mucin_dNdS_less"] = {"U": U, "p_one_sided_CCR_lower": p, "rank_biserial": rbc,
                                       "n_CCR": len(ccr), "n_Mucin": len(mucin)}
    U, p, rbc = mwu(ccr, nonccr, "less")
    tests["CCR_vs_nonCCR_dNdS_less"] = {"U": U, "p_one_sided_CCR_lower": p, "rank_biserial": rbc,
                                        "n_CCR": len(ccr), "n_nonCCR": len(nonccr)}
    U, p, rbc = mwu(cx3c_dnds, mucin, "less")
    tests["CX3C_vs_Mucin_dNdS_less"] = {"U": U, "p_one_sided_CX3C_lower": p, "rank_biserial": rbc,
                                        "n_CX3C": len(cx3c_dnds), "n_Mucin": len(mucin)}

    def counts(reg_codons):
        pur = sum(1 for c in reg_codons if c["beta"] < c["alpha"] and c["p"] <= PCUT)
        return pur, len(reg_codons) - pur
    ccr_codons = [c for c in per_codon if c["region"] == "CCR"]
    mucin_codons = [c for c in per_codon if c["region"] == "Mucin"]
    a_pur, a_not = counts(ccr_codons)
    m_pur, m_not = counts(mucin_codons)
    odds, fisher_p = stats.fisher_exact([[a_pur, a_not], [m_pur, m_not]], alternative="greater")
    tests["fisher_purifying_CCR_gt_Mucin"] = {
        "CCR_purifying": a_pur, "CCR_nonpurifying": a_not,
        "Mucin_purifying": m_pur, "Mucin_nonpurifying": m_not,
        "odds_ratio": float(odds), "p_one_sided_CCR_more_purifying": float(fisher_p),
    }

    out = {
        "label": args.label,
        "regions_locked": fr.RSV_REGIONS_LOCKED,
        "a2_ref_len": ref_len,
        "a2_ref_expected": fr.RSV_G_LENGTH,
        "n_codons_total": len(per_codon),
        "region_definitions": {
            "Nterm": "A2 1-156", "CCR": "A2 157-198", "CX3C_nested": "A2 182-186",
            "Mucin": "A2 199-298", "insertion": "no A2 residue (ON1 duplication)",
        },
        "region_stats": regions,
        "tests": tests,
        "pcut_for_counts": PCUT,
        "metric": "dN - dS  (FEL beta - alpha) per codon",
    }
    outpath = Path("results") / f"{args.label}_region_selection.json"
    json.dump(out, open(outpath, "w"), indent=2)

    csvpath = Path("results") / f"{args.label}_region_selection_percodon.csv"
    with open(csvpath, "w", newline="") as fh:
        fh.write("aln_col,a2,region,alpha,beta,dNdS,p\n")
        for c in per_codon:
            a2v = "" if c["a2"] is None else c["a2"]
            fh.write(f"{c['aln_col']},{a2v},{c['region']},"
                     f"{c['alpha']:.6g},{c['beta']:.6g},{c['dNdS']:.6g},{c['p']:.6g}\n")

    print(f"\n=== REGION-LEVEL SELECTION ({args.label}) ===")
    print(f"A2 ref len {ref_len} (expected {fr.RSV_G_LENGTH}); {len(per_codon)} codons total\n")
    hdr = f"{'region':<12}{'n':>5}{'mean dN-dS':>12}{'median':>10}{'divers.':>9}{'purif.':>8}{'%purif':>8}"
    print(hdr); print("-" * len(hdr))
    order = ["Nterm", "CCR", "CX3C_nested", "Mucin", "insertion"]
    for reg in order:
        s = regions.get(reg)
        if not s:
            continue
        print(f"{reg:<12}{s['n_codons']:>5}{s['mean_dN_minus_dS']:>12.4f}"
              f"{s['median_dN_minus_dS']:>10.4f}{s['n_diversifying_p05']:>9}"
              f"{s['n_purifying_p05']:>8}{s['frac_purifying_p05']*100:>7.1f}%")
    kw = tests.get("kruskal_Nterm_CCR_Mucin")
    if kw:
        print(f"\nOmnibus (Kruskal-Wallis, Nterm/CCR/Mucin): H = {kw['H']:.3f}, "
              f"df = {kw['df']}, p = {kw['p']:.4g}, eps^2 = {kw['epsilon_squared']:.3f}")
    print("\nDirectional tests (Mann-Whitney U on dN-dS; rbc = rank-biserial effect):")
    t = tests["CCR_vs_Mucin_dNdS_less"]
    print(f"  H1 CCR < Mucin (more purifying): p = {t['p_one_sided_CCR_lower']:.4g}"
          f"  rbc = {t['rank_biserial']:+.3f}  (n_CCR={t['n_CCR']}, n_Mucin={t['n_Mucin']})")
    t = tests["CCR_vs_nonCCR_dNdS_less"]
    print(f"     CCR < non-CCR              : p = {t['p_one_sided_CCR_lower']:.4g}"
          f"  rbc = {t['rank_biserial']:+.3f}")
    t = tests["CX3C_vs_Mucin_dNdS_less"]
    print(f"     CX3C < Mucin              : p = {t['p_one_sided_CX3C_lower']:.4g}"
          f"  rbc = {t['rank_biserial']:+.3f}")
    t = tests["fisher_purifying_CCR_gt_Mucin"]
    print(f"  Fisher purifying CCR>Mucin   : OR={t['odds_ratio']:.2f}, "
          f"p = {t['p_one_sided_CCR_more_purifying']:.4g}")
    print(f"\nwrote -> {outpath}")
    print(f"wrote -> {csvpath}")
    return out


if __name__ == "__main__":
    main()
