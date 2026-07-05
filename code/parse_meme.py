import argparse
import csv
import json
import math
from pathlib import Path

import functional_regions as fr
import map_numbering
import map_numbering_rsv


DATASETS = [
    {
        "dataset": "rsv_a",
        "label": "RSV-A G",
        "kind": "rsv",
        "meme_json": "results/meme/rsv_a_meme.json",
        "fel_json": "results/rsv_a_fel.json",
        "coord_aln": "results/rsv_a_prot_aligned_with_A2.fasta",
    },
    {
        "dataset": "rsv_b",
        "label": "RSV-B G",
        "kind": "rsv",
        "meme_json": "results/meme/rsv_b_meme.json",
        "fel_json": "results/rsv_b_fel.json",
        "coord_aln": "results/rsv_b_prot_aligned_with_A2.fasta",
    },
    {
        "dataset": "flu_h3_ha",
        "label": "Influenza H3N2 HA",
        "kind": "flu",
        "meme_json": "results/meme/flu_h3_ha_meme.json",
        "fel_json": "results/flu_h3_ha_fel.json",
        "coord_aln": "results/flu_h3_ha_prot_aligned.fasta",
    },
]


SUMMARY_FIELDS = [
    "dataset",
    "protein",
    "aln_codon_pos",
    "coordinate_system",
    "reference_residue",
    "region",
    "a2_mappable",
    "in_ccr",
    "in_cx3c",
    "in_mucin",
    "fel_positive",
    "fel_alpha_dS",
    "fel_beta_dN",
    "fel_beta_minus_alpha",
    "fel_p_value",
    "meme_positive",
    "meme_alpha",
    "meme_beta_plus",
    "meme_p_plus",
    "meme_omega_plus",
    "meme_lrt",
    "meme_p_value",
    "meme_branches_under_selection",
]


MEME_FIELDS = [
    "dataset",
    "protein",
    "aln_codon_pos",
    "coordinate_system",
    "reference_residue",
    "region",
    "a2_mappable",
    "in_ccr",
    "in_cx3c",
    "in_mucin",
    "alpha",
    "beta_plus",
    "p_plus",
    "omega_plus",
    "lrt",
    "p_value",
    "branches_under_selection",
]


COUNT_FIELDS = [
    "dataset",
    "method",
    "n_positive",
    "n_a2_mappable",
    "n_in_ccr",
    "n_in_cx3c",
    "n_in_mucin",
    "n_on1_or_unmappable",
    "n_other_mapped",
    "n_h3_antigenic",
]


def fmt(value):
    if value in ("", None):
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isinf(value):
            return "inf"
        if math.isnan(value):
            return ""
        return f"{value:.6g}"
    return str(value)


def load_fel_positive(path, pcut):
    data = json.load(open(path, encoding="utf-8"))
    rows = data["MLE"]["content"]["0"]
    out = {}
    for pos, row in enumerate(rows, start=1):
        alpha, beta, p_value = float(row[0]), float(row[1]), float(row[4])
        if p_value <= pcut and beta > alpha:
            out[pos] = {
                "alpha": alpha,
                "beta": beta,
                "beta_minus_alpha": beta - alpha,
                "p_value": p_value,
            }
    return out


def load_meme_positive(path, pcut):
    data = json.load(open(path, encoding="utf-8"))
    rows = data["MLE"]["content"]["0"]
    out = {}
    for pos, row in enumerate(rows, start=1):
        alpha = float(row[0])
        beta_plus = float(row[3])
        p_plus = float(row[4])
        lrt = float(row[5])
        p_value = float(row[6])
        branches = float(row[7])
        if p_value <= pcut and p_plus > 0 and beta_plus > alpha:
            if alpha == 0:
                omega_plus = math.inf if beta_plus > 0 else ""
            else:
                omega_plus = beta_plus / alpha
            out[pos] = {
                "alpha": alpha,
                "beta_plus": beta_plus,
                "p_plus": p_plus,
                "omega_plus": omega_plus,
                "lrt": lrt,
                "p_value": p_value,
                "branches": branches,
            }
    return out


def h3_region(h3_residue):
    if h3_residue is None:
        return "outside mature HA1"
    hits = [
        site for site, residues in sorted(fr.H3_ANTIGENIC_SITES.items())
        if h3_residue in residues
    ]
    if hits:
        return "H3 antigenic site " + "/".join(hits)
    return "non-antigenic HA1"


def rsv_region(a2_residue):
    if a2_residue is None:
        return "ON1/unmappable insertion"
    if a2_residue in fr.RSV_G_CX3C:
        return "CX3C motif"
    if a2_residue in fr.RSV_G_CCR:
        return "CCR"
    if 199 <= a2_residue <= fr.RSV_G_LENGTH:
        return "C-terminal mucin-like domain"
    if 1 <= a2_residue <= 156:
        return "N-terminal region"
    return "outside RSV-A2 G"


def build_mapper(dataset):
    if dataset["kind"] == "rsv":
        aln_to_a2, ref_len = map_numbering_rsv.build_aln_to_a2(
            dataset["coord_aln"], "M11486_A2"
        )

        def mapper(pos):
            a2 = aln_to_a2.get(pos)
            return {
                "coordinate_system": "RSV-A2 mature G",
                "reference_residue": a2 if a2 is not None else "",
                "region": rsv_region(a2),
                "a2_mappable": a2 is not None,
                "in_ccr": a2 in fr.RSV_G_CCR if a2 is not None else False,
                "in_cx3c": a2 in fr.RSV_G_CX3C if a2 is not None else False,
                "in_mucin": 199 <= a2 <= fr.RSV_G_LENGTH if a2 is not None else False,
                "ref_len": ref_len,
            }

        return mapper

    col_to_h3, ref_id = map_numbering.build_column_to_h3(dataset["coord_aln"])

    def mapper(pos):
        h3 = col_to_h3.get(pos)
        return {
            "coordinate_system": "H3 HA1",
            "reference_residue": h3 if h3 is not None else "",
            "region": h3_region(h3),
            "a2_mappable": "",
            "in_ccr": "",
            "in_cx3c": "",
            "in_mucin": "",
            "ref_id": ref_id,
        }

    return mapper


def comparison_rows(dataset, fel, meme):
    mapper = build_mapper(dataset)
    rows = []
    for pos in sorted(set(fel) | set(meme)):
        meta = mapper(pos)
        f = fel.get(pos, {})
        m = meme.get(pos, {})
        rows.append({
            "dataset": dataset["dataset"],
            "protein": dataset["label"],
            "aln_codon_pos": pos,
            "coordinate_system": meta["coordinate_system"],
            "reference_residue": meta["reference_residue"],
            "region": meta["region"],
            "a2_mappable": meta["a2_mappable"],
            "in_ccr": meta["in_ccr"],
            "in_cx3c": meta["in_cx3c"],
            "in_mucin": meta["in_mucin"],
            "fel_positive": pos in fel,
            "fel_alpha_dS": f.get("alpha", ""),
            "fel_beta_dN": f.get("beta", ""),
            "fel_beta_minus_alpha": f.get("beta_minus_alpha", ""),
            "fel_p_value": f.get("p_value", ""),
            "meme_positive": pos in meme,
            "meme_alpha": m.get("alpha", ""),
            "meme_beta_plus": m.get("beta_plus", ""),
            "meme_p_plus": m.get("p_plus", ""),
            "meme_omega_plus": m.get("omega_plus", ""),
            "meme_lrt": m.get("lrt", ""),
            "meme_p_value": m.get("p_value", ""),
            "meme_branches_under_selection": m.get("branches", ""),
        })
    return rows


def meme_rows(dataset, meme):
    mapper = build_mapper(dataset)
    rows = []
    for pos in sorted(meme):
        meta = mapper(pos)
        m = meme[pos]
        rows.append({
            "dataset": dataset["dataset"],
            "protein": dataset["label"],
            "aln_codon_pos": pos,
            "coordinate_system": meta["coordinate_system"],
            "reference_residue": meta["reference_residue"],
            "region": meta["region"],
            "a2_mappable": meta["a2_mappable"],
            "in_ccr": meta["in_ccr"],
            "in_cx3c": meta["in_cx3c"],
            "in_mucin": meta["in_mucin"],
            "alpha": m["alpha"],
            "beta_plus": m["beta_plus"],
            "p_plus": m["p_plus"],
            "omega_plus": m["omega_plus"],
            "lrt": m["lrt"],
            "p_value": m["p_value"],
            "branches_under_selection": m["branches"],
        })
    return rows


def count_row(dataset, method, positive_rows):
    if dataset["kind"] == "rsv":
        n_positive = len(positive_rows)
        n_a2 = sum(r["a2_mappable"] is True for r in positive_rows)
        n_ccr = sum(r["in_ccr"] is True for r in positive_rows)
        n_cx3c = sum(r["in_cx3c"] is True for r in positive_rows)
        n_mucin = sum(r["in_mucin"] is True for r in positive_rows)
        n_unmapped = sum(r["a2_mappable"] is False for r in positive_rows)
        return {
            "dataset": dataset["dataset"],
            "method": method,
            "n_positive": n_positive,
            "n_a2_mappable": n_a2,
            "n_in_ccr": n_ccr,
            "n_in_cx3c": n_cx3c,
            "n_in_mucin": n_mucin,
            "n_on1_or_unmappable": n_unmapped,
            "n_other_mapped": n_a2 - n_ccr - n_mucin,
            "n_h3_antigenic": "",
        }

    n_antigenic = sum("H3 antigenic site" in r["region"] for r in positive_rows)
    return {
        "dataset": dataset["dataset"],
        "method": method,
        "n_positive": len(positive_rows),
        "n_a2_mappable": "",
        "n_in_ccr": "",
        "n_in_cx3c": "",
        "n_in_mucin": "",
        "n_on1_or_unmappable": "",
        "n_other_mapped": "",
        "n_h3_antigenic": n_antigenic,
    }


def write_csv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def markdown_table(rows, fields):
    lines = []
    lines.append("| " + " | ".join(fields) + " |")
    lines.append("| " + " | ".join(["---"] * len(fields)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(field, "")) for field in fields) + " |")
    return "\n".join(lines)


def write_markdown(path, comparison, counts, pcut):
    detail_fields = [
        "dataset",
        "aln_codon_pos",
        "reference_residue",
        "region",
        "a2_mappable",
        "fel_positive",
        "fel_beta_dN",
        "fel_p_value",
        "meme_positive",
        "meme_beta_plus",
        "meme_p_plus",
        "meme_omega_plus",
        "meme_p_value",
    ]
    count_fields = [
        "dataset",
        "method",
        "n_positive",
        "n_a2_mappable",
        "n_in_ccr",
        "n_in_cx3c",
        "n_in_mucin",
        "n_on1_or_unmappable",
        "n_other_mapped",
        "n_h3_antigenic",
    ]
    text = [
        f"# MEME episodic-selection sensitivity summary (p <= {pcut:g})",
        "",
        "Raw MEME JSON files are in `results/meme/*_meme.json`; this summary is generated by `code/parse_meme.py`.",
        "",
        "## Positive-site counts",
        "",
        markdown_table(counts, count_fields),
        "",
        "## FEL/MEME positive-site comparison",
        "",
        markdown_table(comparison, detail_fields),
        "",
    ]
    path.write_text("\n".join(text), encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--out-dir", default="results/meme")
    args = parser.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_comparison = []
    all_counts = []
    for dataset in DATASETS:
        fel = load_fel_positive(dataset["fel_json"], args.p_cutoff)
        meme = load_meme_positive(dataset["meme_json"], args.p_cutoff)
        cmp_rows = comparison_rows(dataset, fel, meme)
        meme_site_rows = meme_rows(dataset, meme)
        fel_site_rows = [row for row in cmp_rows if row["fel_positive"]]

        write_csv(
            out_dir / f"{dataset['dataset']}_meme_sites.csv",
            meme_site_rows,
            MEME_FIELDS,
        )

        all_comparison.extend(cmp_rows)
        all_counts.append(count_row(dataset, "FEL", fel_site_rows))
        all_counts.append(count_row(dataset, "MEME", meme_site_rows))

    write_csv(out_dir / "summary_fel_meme_sites.csv", all_comparison, SUMMARY_FIELDS)
    write_csv(out_dir / "summary_meme_region_counts.csv", all_counts, COUNT_FIELDS)
    write_markdown(
        out_dir / "summary_fel_meme_sites.md",
        all_comparison,
        all_counts,
        args.p_cutoff,
    )

    print(f"p-value cutoff: {args.p_cutoff:g}")
    print(f"wrote -> {out_dir / 'summary_fel_meme_sites.csv'}")
    print(f"wrote -> {out_dir / 'summary_meme_region_counts.csv'}")
    print(f"wrote -> {out_dir / 'summary_fel_meme_sites.md'}")


if __name__ == "__main__":
    main()
