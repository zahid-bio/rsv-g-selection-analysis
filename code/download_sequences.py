import sys
import time
import random
from pathlib import Path

from Bio import Entrez, SeqIO

CONFIG = {
    "name": "flu_h3_ha_nt",
    "query": (
        '"Influenza A virus"[Organism] AND H3N2[All Fields] '
        'AND hemagglutinin[Title] AND "complete cds"[Title] '
        'AND "Homo sapiens"[All Fields]'
    ),
    "year_bins": [
        (1968, 1998), (1998, 2004), (2004, 2008), (2008, 2012),
        (2012, 2016), (2016, 2019), (2019, 2022), (2022, 2024),
        (2024, 2027),
    ],
    "per_bin": 22,
    "pool_per_bin": 400,
    "min_len": 1650,
    "max_len": 1780,
    "random_seed": 42,
}

PROJECT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT / "data"
NOTES_DIR = PROJECT / "notes"
SECRETS = PROJECT / "code" / ".entrez_secrets.txt"


def load_secrets():
    if not SECRETS.exists():
        sys.exit(f"ERROR: secrets file not found at {SECRETS}")
    secrets = {}
    for line in SECRETS.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        secrets[k.strip()] = v.strip()
    if not secrets.get("NCBI_EMAIL"):
        sys.exit("ERROR: NCBI_EMAIL empty in code/.entrez_secrets.txt. Add it and re-run.")
    return secrets.get("NCBI_EMAIL"), secrets.get("NCBI_API_KEY", "")


def esearch_window(query, y0, y1, retmax):
    h = Entrez.esearch(
        db="nuccore", term=query, retmax=retmax,
        datetype="pdat", mindate=f"{y0}/01/01", maxdate=f"{y1-1}/12/31",
    )
    res = Entrez.read(h); h.close()
    return res["IdList"], int(res["Count"])


def efetch_fasta(ids):
    recs = []
    for start in range(0, len(ids), 100):
        batch = ids[start:start + 100]
        for attempt in range(3):
            try:
                h = Entrez.efetch(db="nuccore", id=",".join(batch),
                                  rettype="fasta", retmode="text")
                recs.extend(list(SeqIO.parse(h, "fasta")))
                h.close()
                break
            except Exception as e:
                print(f"      batch {start} attempt {attempt+1} failed: {e}")
                time.sleep(3)
        time.sleep(0.4)
    return recs


def main():
    email, api_key = load_secrets()
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    random.seed(CONFIG["random_seed"])
    name = CONFIG["name"]

    print(f"[1/4] Stratified search across {len(CONFIG['year_bins'])} date windows")
    chosen_ids = []
    for (y0, y1) in CONFIG["year_bins"]:
        ids, count = esearch_window(CONFIG["query"], y0, y1, CONFIG["pool_per_bin"])
        take = min(CONFIG["per_bin"], len(ids))
        picked = random.sample(ids, take) if take else []
        chosen_ids.extend(picked)
        print(f"      {y0}-{y1-1}: {count} matches, pool {len(ids)}, took {len(picked)}")
        time.sleep(0.3)
    chosen_ids = list(dict.fromkeys(chosen_ids))
    print(f"      total IDs selected: {len(chosen_ids)}")

    print(f"[2/4] Fetching {len(chosen_ids)} FASTA records ...")
    records = efetch_fasta(chosen_ids)
    print(f"      fetched {len(records)} raw records.")

    kept, seen = [], set()
    for rec in records:
        L = len(rec.seq)
        if not (CONFIG["min_len"] <= L <= CONFIG["max_len"]):
            continue
        s = str(rec.seq).upper()
        if s in seen:
            continue
        seen.add(s); kept.append(rec)
    print(f"[3/4] After length QC + dedupe: {len(kept)} unique records.")

    DATA_DIR.mkdir(exist_ok=True); NOTES_DIR.mkdir(exist_ok=True)
    fasta_out = DATA_DIR / f"{name}.fasta"
    tsv_out = NOTES_DIR / f"{name}_accessions.tsv"
    SeqIO.write(kept, fasta_out, "fasta")
    with open(tsv_out, "w", encoding="utf-8") as f:
        f.write("accession\tlength_nt\tdescription\n")
        for rec in kept:
            f.write(f"{rec.id}\t{len(rec.seq)}\t{rec.description}\n")
    print(f"[4/4] Wrote {len(kept)} sequences -> {fasta_out}")
    print(f"      metadata -> {tsv_out}")
    print("Next: extract_cds.py -> codon_align.py -> IQ-TREE -> FEL.")


if __name__ == "__main__":
    main()
