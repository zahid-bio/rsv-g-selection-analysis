import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def longest_orf(seq, min_nt, max_nt):
    best = None
    for strand, nuc in [("+", seq), ("-", seq.reverse_complement())]:
        s = str(nuc).upper()
        for frame in range(3):
            i = frame
            while i < len(s) - 2:
                if s[i:i+3] == "ATG":
                    j = i
                    while j < len(s) - 2:
                        codon = s[j:j+3]
                        if codon in ("TAA", "TAG", "TGA"):
                            orf = s[i:j]
                            if min_nt <= len(orf) <= max_nt and len(orf) % 3 == 0:
                                if best is None or len(orf) > len(best):
                                    best = orf
                            break
                        j += 3
                    i += 3
                else:
                    i += 3
    return Seq(best) if best else None


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    inp, outp = Path(sys.argv[1]), Path(sys.argv[2])
    min_codons = int(sys.argv[3]) if len(sys.argv) > 3 else 540
    max_codons = int(sys.argv[4]) if len(sys.argv) > 4 else 590
    min_nt, max_nt = min_codons * 3, max_codons * 3

    recs = list(SeqIO.parse(inp, "fasta"))
    if not recs:
        sys.exit("Error: Input FASTA file is empty or does not contain valid records.")
    kept, dropped = [], []
    for r in recs:
        orf = longest_orf(r.seq, min_nt, max_nt)
        if orf is None:
            dropped.append(r.id)
            continue
        prot = str(orf.translate())
        if "*" in prot:
            dropped.append(r.id)
            continue
        kept.append(SeqRecord(orf, id=r.id, description=r.description))

    SeqIO.write(kept, outp, "fasta")
    print(f"input records      : {len(recs)}")
    print(f"clean CDS extracted: {len(kept)}  -> {outp}")
    print(f"dropped (no ORF)   : {len(dropped)}")
    if dropped:
        print("  dropped ids:", ", ".join(dropped[:15]) + (" ..." if len(dropped) > 15 else ""))
    if kept:
        lens = [len(r.seq) for r in kept]
        print(f"CDS length range   : {min(lens)}-{max(lens)} nt "
              f"({min(lens)//3}-{max(lens)//3} codons); all divisible by 3: "
              f"{all(L % 3 == 0 for L in lens)}")


if __name__ == "__main__":
    main()
