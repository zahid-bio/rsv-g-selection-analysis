# Reproducibility Note

## Frozen Inputs

- RSV-A input: `data/rsv_a_nt.fasta`, 108 records.
- RSV-B input: `data/rsv_b_nt.fasta`, 115 records.
- Accession lists: `notes/rsv_a_nt_accessions.tsv` and `notes/rsv_b_nt_accessions.tsv`.
- The exact original RSV Entrez query text and random seed were not captured at download time.
- Therefore, RSV sequence acquisition is documented as a frozen-input limitation. The analysis is reproducible from the archived FASTA/accession files forward, but not as a byte-for-byte re-download.

## Region Definitions

- RSV-A2 G reference: GenBank M11486, 298 aa.
- CCR: RSV-A2 residues 157-198.
- CX3C motif: RSV-A2 residues 182-186.
- C-terminal mucin-like domain: RSV-A2 residues 199-298.
- The region definitions were specified and committed before RSV selected-site positions were inspected. This should be described as "pre-specified," not as formal preregistration unless an immutable external preregistration artifact is added.

## Main Analysis Files

- FEL raw JSON: `results/rsv_a_fel.json`, `results/rsv_b_fel.json`, `results/flu_h3_ha_fel.json`.
- FEL parsed summaries: `results/*_fel_sites.csv`, `results/summary_selected_sites.csv`.
- Enrichment summaries: `results/summary_enrichment.md`, `results/summary_enrichment.csv`.
- Region-level summaries: `results/summary_region_selection.csv`, `results/*_region_selection.json`, `results/*_region_selection_percodon.csv`.

## MEME Sensitivity Analysis

### Raw Outputs

Raw MEME outputs were generated using HyPhy in Docker and are stored under `results/meme/`.

**HyPhy version:** 2.5.100 (Docker image: `quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0`)
**Analysis version:** MEME 4.1 as reported in console logs
**Genetic code:** Universal
**P-value threshold:** 0.05 (unadjusted, screening-level)
**Branches included:** All branches

### Commands to Regenerate MEME JSON Files

Run from repository root (where `results/` is a subdirectory):

```powershell
# RSV-A G
docker run --rm --mount type=bind,source="${PWD}",target=/data quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0 hyphy meme --alignment /data/results/rsv_a_codon_aligned.fasta --tree /data/results/rsv_a_codon_iqtree.treefile --code Universal --pvalue 0.05 --output /data/results/meme/rsv_a_meme.json

# RSV-B G
docker run --rm --mount type=bind,source="${PWD}",target=/data quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0 hyphy meme --alignment /data/results/rsv_b_codon_aligned.fasta --tree /data/results/rsv_b_codon_iqtree.treefile --code Universal --pvalue 0.05 --output /data/results/meme/rsv_b_meme.json

# Influenza H3N2 HA (positive control)
docker run --rm --mount type=bind,source="${PWD}",target=/data quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0 hyphy meme --alignment /data/results/flu_h3_ha_codon_aligned.fasta --tree /data/results/flu_h3_ha_codon_iqtree.treefile --code Universal --pvalue 0.05 --output /data/results/meme/flu_h3_ha_meme.json
```

**Input alignments:** Codon-aware alignments created from CDS files with MAFFT and back-translated to codon coordinates.
**Input trees:** Maximum-likelihood trees inferred with IQ-TREE v3.1.3, unrooted, used as-is.
**Output format:** JSON (native HyPhy MLE output).

### Parsed MEME Outputs

The raw JSON files are parsed into tabular format by `code/parse_meme.py`:

```powershell
python code/parse_meme.py --p-cutoff 0.05 --out-dir results/meme
```

**Parsed output files:**
- `results/meme/rsv_a_meme_sites.csv` — MEME-positive sites in RSV-A G
- `results/meme/rsv_b_meme_sites.csv` — MEME-positive sites in RSV-B G
- `results/meme/flu_h3_ha_meme_sites.csv` — MEME-positive sites in influenza HA (positive control)
- `results/meme/summary_meme_region_counts.csv` — Count summary by region (Table S2, top panel)
- `results/meme/summary_fel_meme_sites.csv` — Combined FEL/MEME comparison (Table S2, detailed)
- `results/meme/summary_fel_meme_sites.md` — Markdown-formatted summary

### SHA256 Checksums

All raw JSON and parsed CSV files are checksummed in `results/meme/SHA256SUMS`:

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

Verify with: `cd results/meme && sha256sum -c SHA256SUMS`

### Interpretation Guardrails

- **FEL remains the primary analysis:** Pervasive/site-wise diversifying selection is reported via FEL.
- **MEME is secondary and episodic:** MEME detects episodic diversifying selection on a subset of branches.
- **CX3C/cystine-noose constraint:** MEME did not identify episodic diversifying selection within the CX3C motif itself in this dataset.
- **Broader CCR sensitivity:** MEME did identify one non-CX3C broader-CCR site in each RSV subtype (RSV-A residue 178, RSV-B residue 164); the broader CCR should not be described as completely spared from episodic signal.
- **Not FEL-positive:** The MEME-only CCR sites are not significant in FEL; they represent branch-sparse episodic signals only.

## Statistical Interpretation Notes

- FEL and MEME site-wise p-values use unadjusted p <= 0.05 screening thresholds.
- The selected-site enrichment test is one-sided for enrichment. A result of p = 1.0 means "not enriched," not formal depletion and not proof of absence.
- The pooled RSV-A+B selected-site table is descriptive because RSV-A and RSV-B were analyzed with subtype-specific alignments and trees.
- The region-level dN-dS analysis is threshold-free and uses all mapped codons, but it treats per-codon FEL point estimates as comparable summaries and does not propagate uncertainty in individual estimates or tree/model uncertainty.

## Multiple-testing correction (FDR)

- Benjamini-Hochberg FDR q-values are computed across all codons per dataset and method by `code/apply_fdr.py`, writing `results/summary_fdr_qvalues.csv`. FDR is controlled at q <= 0.05.
- Result: after FDR correction, **no RSV G site remains significant** (FEL or MEME) in either subtype. Only 3 influenza control sites survive: FEL codons 161 and 287 (antigenic sites A/E) and MEME codon 157 (antigenic site A).
- Reproduce: `python code/apply_fdr.py --p-cutoff 0.05 --q-cutoff 0.05 --out results/summary_fdr_qvalues.csv`.
- Interpretation: RSV site-level calls are screening-level; the study's inference rests on the threshold-free region-level analysis (unaffected by any per-site cutoff) and on the influenza control, whose antigenic signal survives FDR.

## Recombination screen (GARD)

- Each RSV G codon alignment was screened with HyPhy GARD (`code/run_gard.sh`).
- Command (MPI; the `TOLERATE_NUMERICAL_ERRORS` environment flag is required for numerical stability with these alignments, otherwise GARD aborts with a `ComputeBranchCache` internal error):
  ```
  docker run --rm --shm-size=1g -v "$PWD":/data quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0 \
    mpirun --allow-run-as-root --mca btl tcp,self --oversubscribe -np 4 HYPHYMPI gard \
    --alignment /data/results/rsv_a_codon_aligned.fasta --output /data/results/rsv_a_gard.json \
    ENV="TOLERATE_NUMERICAL_ERRORS=1;"
  ```
- Result: the exhaustive single-breakpoint scan found **no breakpoint improving the c-AIC** over the no-breakpoint model in either RSV-A (373 candidate breakpoints) or RSV-B (360 candidate breakpoints), and the multi-breakpoint genetic-algorithm search returned no c-AIC improvement. No evidence of recombination; single-tree FEL/MEME analyses are appropriate. Console logs: `results/rsv_a_gard_console.log`, `results/rsv_b_gard_console.log`.
