import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def translate(cds_path, prot_path):
    prot_recs = []
    for r in SeqIO.parse(cds_path, "fasta"):
        if len(r.seq) % 3 != 0:
            raise ValueError(f"Sequence {r.id} length is not a multiple of 3 (length={len(r.seq)}).")
        prot = r.seq.translate()
        prot_recs.append(SeqRecord(prot, id=r.id, description=""))
    SeqIO.write(prot_recs, prot_path, "fasta")
    print(f"translated {len(prot_recs)} CDS -> {prot_path}")


def backtranslate(cds_path, aln_prot_path, out_path):
    cds = {r.id: str(r.seq).upper() for r in SeqIO.parse(cds_path, "fasta")}
    out = []
    for pr in SeqIO.parse(aln_prot_path, "fasta"):
        nt = cds[pr.id]
        aa_len = len(str(pr.seq).replace("-", ""))
        if len(nt) // 3 != aa_len:
            raise ValueError(f"Aligned protein length without gaps ({aa_len}) does not match nucleotide sequence length / 3 ({len(nt) // 3}) for sequence {pr.id}")
        codons, k = [], 0
        for aa in str(pr.seq):
            if aa == "-":
                codons.append("---")
            else:
                if k + 3 > len(nt):
                    raise ValueError(f"Codon index {k} is out of bounds for sequence {pr.id} (length {len(nt)})")
                codons.append(nt[k:k+3])
                k += 3
        out.append(SeqRecord(Seq("".join(codons)), id=pr.id, description=""))
    SeqIO.write(out, out_path, "fasta")
    L = len(out[0].seq)
    print(f"back-translated {len(out)} seqs -> {out_path}")
    print(f"codon alignment length: {L} nt  (divisible by 3: {L % 3 == 0})")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    cmd = sys.argv[1]
    if cmd == "translate" and len(sys.argv) == 4:
        translate(sys.argv[2], sys.argv[3])
    elif cmd == "backtranslate" and len(sys.argv) == 5:
        backtranslate(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
