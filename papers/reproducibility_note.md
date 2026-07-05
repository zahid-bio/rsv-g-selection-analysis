# Reproducibility Note

This note records the reproducibility status of the RSV G selection analysis after
the critical review of `main (5).pdf`.

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

Raw MEME outputs were added in commit `61f475d` and are stored under `results/meme/`.

Commands used from the repository root:

```powershell
docker run --rm --mount type=bind,source="${PWD}",target=/data quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0 hyphy meme --alignment /data/results/rsv_a_codon_aligned.fasta --tree /data/results/rsv_a_codon_iqtree.treefile --code Universal --pvalue 0.05 --output /data/results/meme/rsv_a_meme.json
docker run --rm --mount type=bind,source="${PWD}",target=/data quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0 hyphy meme --alignment /data/results/rsv_b_codon_aligned.fasta --tree /data/results/rsv_b_codon_iqtree.treefile --code Universal --pvalue 0.05 --output /data/results/meme/rsv_b_meme.json
docker run --rm --mount type=bind,source="${PWD}",target=/data quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0 hyphy meme --alignment /data/results/flu_h3_ha_codon_aligned.fasta --tree /data/results/flu_h3_ha_codon_iqtree.treefile --code Universal --pvalue 0.05 --output /data/results/meme/flu_h3_ha_meme.json
python code/parse_meme.py --p-cutoff 0.05 --out-dir results/meme
```

Parsed MEME outputs:

- `results/meme/rsv_a_meme_sites.csv`
- `results/meme/rsv_b_meme_sites.csv`
- `results/meme/flu_h3_ha_meme_sites.csv`
- `results/meme/summary_meme_region_counts.csv`
- `results/meme/summary_fel_meme_sites.csv`
- `results/meme/summary_fel_meme_sites.md`

Interpretation guardrail:

- FEL remains the primary pervasive/site-wise analysis.
- MEME is a secondary episodic/site-wise sensitivity analysis.
- MEME did not identify episodic diversifying selection in the CX3C motif in this dataset.
- MEME did identify one non-CX3C broader-CCR site in each RSV subtype, so the broader CCR should not be described as completely spared.

## Statistical Interpretation Notes

- FEL and MEME site-wise p-values use unadjusted p <= 0.05 screening thresholds.
- The selected-site enrichment test is one-sided for enrichment. A result of p = 1.0 means "not enriched," not formal depletion and not proof of absence.
- The pooled RSV-A+B selected-site table is descriptive because RSV-A and RSV-B were analyzed with subtype-specific alignments and trees.
- The region-level dN-dS analysis is threshold-free and uses all mapped codons, but it treats per-codon FEL point estimates as comparable summaries and does not propagate uncertainty in individual estimates or tree/model uncertainty.
- No recombination screen or alternative-tree sensitivity analysis has been run.
