# MEME Analysis Reproducibility — Quick Start

**Status:** ✅ Fully reproducible  
**Location:** `results/meme/`  
**Last verified:** 2026-07-04

---

## Verify in 3 Steps

```powershell
# Step 1: Check raw MEME JSON files are valid
python -c "
import json
for f in ['rsv_a', 'rsv_b', 'flu_h3_ha']:
    json.load(open(f'results/meme/{f}_meme.json'))
    print(f'{f}: ✅')
"

# Step 2: Regenerate CSV summaries from JSON
python code/parse_meme.py --p-cutoff 0.05 --out-dir results/meme

# Step 3: Verify checksums match
cd results/meme
sha256sum -c SHA256SUMS
```

If all three steps pass with "OK" messages, reproducibility is confirmed.

---

## What You Get

### Raw Outputs (Archived, Non-regenerable)
- `rsv_a_meme.json` (1.2 MB) — HyPhy MEME output for RSV-A G
- `rsv_b_meme.json` (1.3 MB) — HyPhy MEME output for RSV-B G
- `flu_h3_ha_meme.json` (2.2 MB) — HyPhy MEME output for influenza HA

### Regenerable Outputs (CSV summaries)
- `summary_meme_region_counts.csv` — Count by region (Table S2 panel A)
- `summary_fel_meme_sites.csv` — FEL/MEME comparison (Table S2 panel B)
- `rsv_a_meme_sites.csv` — RSV-A MEME-positive sites with parameters
- `rsv_b_meme_sites.csv` — RSV-B MEME-positive sites with parameters
- `flu_h3_ha_meme_sites.csv` — Influenza MEME-positive sites

---

## Key Results (Validation)

| Claim | Result | Status |
|-------|--------|--------|
| CX3C motif is constrained | 0 MEME sites in CX3C (RSV-A & B) | ✅ Match manuscript |
| Broader CCR has sparse episodic signal | 1 site per subtype (non-CX3C CCR) | ✅ Match manuscript |
| Mucin domain is hotspot | 4-5 MEME sites per subtype | ✅ Match manuscript |
| Influenza validates pipeline | 6/11 sites in antigenic regions | ✅ Match manuscript |

---

## HyPhy Specification

- **Version:** 2.5.100 (Docker: `quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0`)
- **Genetic code:** Universal
- **P-value threshold:** 0.05 (unadjusted, screening-level)
- **Branches:** All (no subsetting)

---

## Commands to Re-run MEME (if needed)

From repository root:

```bash
# RSV-A G
docker run --rm --mount type=bind,source="${PWD}",target=/data \
  quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0 hyphy meme \
  --alignment /data/results/rsv_a_codon_aligned.fasta \
  --tree /data/results/rsv_a_codon_iqtree.treefile \
  --code Universal --pvalue 0.05 \
  --output /data/results/meme/rsv_a_meme.json

# RSV-B G
docker run --rm --mount type=bind,source="${PWD}",target=/data \
  quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0 hyphy meme \
  --alignment /data/results/rsv_b_codon_aligned.fasta \
  --tree /data/results/rsv_b_codon_iqtree.treefile \
  --code Universal --pvalue 0.05 \
  --output /data/results/meme/rsv_b_meme.json

# Influenza H3N2 HA
docker run --rm --mount type=bind,source="${PWD}",target=/data \
  quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0 hyphy meme \
  --alignment /data/results/flu_h3_ha_codon_aligned.fasta \
  --tree /data/results/flu_h3_ha_codon_iqtree.treefile \
  --code Universal --pvalue 0.05 \
  --output /data/results/meme/flu_h3_ha_meme.json
```

Then regenerate CSVs: `python code/parse_meme.py --p-cutoff 0.05 --out-dir results/meme`

---

## SHA256 Checksums (to verify)

```
0E3B30E989C9F5FBA0E4B1BBF30CE77FD2DCAAC2F12E786B99C7ABAEC15407D6  flu_h3_ha_meme.json
5B7F524A641FC6DCA2F4164794E4930D62B1E3F6D9FBE5B17AF1B625D7CA60D6  rsv_a_meme.json
2990A19CE7ADBBD8168F7EB94AF4AC067BE58FDBC35D721D2331385AB153ED0B  rsv_b_meme.json
77C9F5E39F79B0B5378DDE633F646301F95030EFCD531AF3FA2D54D30E36AE42  flu_h3_ha_meme_sites.csv
D80B417105B082F071C1A4FBD7ECE078281C858AD707D11893AA8A0ACE1A8D64  rsv_a_meme_sites.csv
3B4E453BAD1AE04EB3FDD55504179507040DED29767F8B738245B2502B32DB17  rsv_b_meme_sites.csv
333E87B7DD85D6FA25761BC47E516B3839587A8D22A0F796DE2FC43AF4441BB0  summary_fel_meme_sites.csv
CC69757E0A1CB1BB5B22CC9970C6B6041E6E86E4E78AF5343D293806F9F855B5  summary_meme_region_counts.csv
```

---

## Documentation

- **Full details:** `papers/reproducibility_note.md`
- **Verification steps:** `notes/MEME_reproducibility_checklist.md`
- **Comprehensive guide:** `notes/MEME_manuscript_reproducibility_summary.md`
- **This document:** `MEME_REPRODUCIBILITY_QUICKSTART.md`

---

## Questions?

1. **Why are my regenerated CSVs different sizes?** — They should be identical. Check that `code/parse_meme.py` ran to completion without errors.

2. **Can I re-run MEME from the Docker command?** — Yes, if you have Docker installed. Use the commands above. Results should match checksums within floating-point precision.

3. **What if MEME gives different p-values?** — Minor differences (p < 0.0001) may occur due to HyPhy numerical precision or random seed. Use checksummed JSON files as ground truth.

4. **Are the results in the manuscript correct?** — Yes. CSV counts and site locations match Table S2. No changes needed.

---

*Reproducibility verified 2026-07-04*
