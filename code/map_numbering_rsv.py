import sys, csv, argparse
from pathlib import Path
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
import functional_regions as fr

CCR = set(fr.RSV_G_CCR)
CX3C = set(fr.RSV_G_CX3C)


def build_aln_to_a2(aln_path, ref_id):
    ref = None
    for r in SeqIO.parse(aln_path, "fasta"):
        if r.id == ref_id:
            ref = r
            break
    if ref is None:
        raise SystemExit(f"reference id {ref_id!r} not found in {aln_path}")
    mapping = {}
    a2 = 0
    for col, aa in enumerate(str(ref.seq), start=1):
        if aa == "-":
            mapping[col] = None
        else:
            a2 += 1
            mapping[col] = a2
    total_ref_len = a2
    if total_ref_len != fr.RSV_G_LENGTH:
        print(f"WARN: A2 reference non-gap length is {total_ref_len}, "
              f"but functional_regions.RSV_G_LENGTH is {fr.RSV_G_LENGTH}", file=sys.stderr)
    return mapping, total_ref_len


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln-with-A2", required=True)
    ap.add_argument("--fel-sites-csv", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--ref-id", default="M11486_A2")
    args = ap.parse_args(argv)

    aln_to_a2, ref_len = build_aln_to_a2(args.aln_with_A2, args.ref_id)

    rows = list(csv.DictReader(open(args.fel_sites_csv)))
    mapped = []
    for r in rows:
        col = int(r["aln_codon_pos"])
        a2 = aln_to_a2.get(col)
        mapped.append({
            "aln_codon_pos": col,
            "a2_residue": a2 if a2 is not None else "",
            "in_ccr": (a2 in CCR) if a2 is not None else "",
            "in_cx3c": (a2 in CX3C) if a2 is not None else "",
            "beta_dN": r.get("beta_dN", ""),
            "p_value": r.get("p_value", ""),
        })

    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mapped[0].keys()))
        w.writeheader()
        for m in mapped:
            w.writerow(m)

    n_total = len(mapped)
    n_a2 = sum(1 for m in mapped if m["a2_residue"] != "")
    n_outside = n_total - n_a2
    n_ccr = sum(1 for m in mapped if m["in_ccr"] is True)
    n_cx3c = sum(1 for m in mapped if m["in_cx3c"] is True)

    print(f"reference id            : {args.ref_id}")
    print(f"A2 reference non-gap len : {ref_len}  (expected {fr.RSV_G_LENGTH})")
    print(f"input FEL rows          : {n_total}")
    print(f"  mapped to A2 residue   : {n_a2}")
    print(f"  in ON1-insertion (no A2): {n_outside}")
    print(f"  in CCR ({min(CCR)}-{max(CCR)}) : {n_ccr}")
    print(f"  in CX3C ({min(CX3C)}-{max(CX3C)}): {n_cx3c}")
    a2_nums = sorted(int(m["a2_residue"]) for m in mapped if m["a2_residue"] != "")
    print(f"  A2 residue numbers     : {a2_nums}")
    print(f"wrote -> {args.out_csv}")


if __name__ == "__main__":
    main()
