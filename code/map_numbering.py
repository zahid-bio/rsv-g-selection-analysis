import sys, csv
from Bio import SeqIO

SIGNAL_PEPTIDE = 16
HA1_LENGTH = 329


def build_column_to_h3(prot_aln_path):
    recs = list(SeqIO.parse(prot_aln_path, "fasta"))
    ref = min(recs, key=lambda r: str(r.seq).count("-"))
    col_to_h3 = {}
    precursor_idx = 0
    for col, aa in enumerate(str(ref.seq), start=1):
        if aa == "-":
            col_to_h3[col] = None
            continue
        precursor_idx += 1
        h3 = precursor_idx - SIGNAL_PEPTIDE
        col_to_h3[col] = h3 if 1 <= h3 <= HA1_LENGTH else None
    return col_to_h3, ref.id


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    prot_aln, fel_csv, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]
    col_to_h3, ref_id = build_column_to_h3(prot_aln)

    rows = list(csv.DictReader(open(fel_csv)))
    mapped, in_ha1, outside = [], 0, 0
    for r in rows:
        col = int(r["aln_codon_pos"])
        h3 = col_to_h3.get(col)
        if h3 is None:
            outside += 1
        else:
            in_ha1 += 1
        mapped.append({"aln_codon_pos": col, "h3_ha1_number": h3,
                       "beta_dN": r["beta_dN"], "p_value": r["p_value"]})

    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["aln_codon_pos", "h3_ha1_number",
                                          "beta_dN", "p_value"])
        w.writeheader()
        for m in mapped:
            w.writerow(m)

    print(f"reference sequence used : {ref_id} (fewest gaps)")
    print(f"signal peptide offset   : {SIGNAL_PEPTIDE}  (mature HA1 = precursor - 16)")
    print(f"selected sites total    : {len(rows)}")
    print(f"  mapped into HA1        : {in_ha1}")
    print(f"  outside HA1 (SP/HA2)   : {outside}")
    h3nums = sorted(m["h3_ha1_number"] for m in mapped if m["h3_ha1_number"])
    print(f"  H3 HA1 residue numbers : {h3nums}")
    print(f"wrote -> {out_csv}")


if __name__ == "__main__":
    main()
