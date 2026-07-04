import sys, json, csv

def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    fel_json, out_csv = sys.argv[1], sys.argv[2]
    pcut = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05

    d = json.load(open(fel_json))
    rows = d["MLE"]["content"]["0"]

    pos, neg = [], []
    for i, r in enumerate(rows, start=1):
        alpha, beta, _, lrt, p, bl = r[0], r[1], r[2], r[3], r[4], r[5]
        if p <= pcut and beta > alpha:
            pos.append((i, alpha, beta, beta - alpha, p))
        elif p <= pcut and beta < alpha:
            neg.append(i)

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["aln_codon_pos", "alpha_dS", "beta_dN", "dN_minus_dS", "p_value"])
        for row in pos:
            w.writerow([row[0], f"{row[1]:.4f}", f"{row[2]:.4f}",
                        f"{row[3]:.4f}", f"{row[4]:.5f}"])

    print(f"FEL file        : {fel_json}")
    print(f"total codon sites: {len(rows)}")
    print(f"p-value cutoff  : {pcut}")
    print(f"POSITIVE (diversifying) sites: {len(pos)}")
    print(f"PURIFYING sites              : {len(neg)}")
    print(f"positively-selected alignment positions:")
    print("  " + ", ".join(str(p[0]) for p in pos))
    print(f"wrote -> {out_csv}")

if __name__ == "__main__":
    main()
