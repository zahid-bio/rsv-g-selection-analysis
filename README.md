# Selection analysis of the RSV G glycoprotein central conserved region

Reproducible pipeline and data for the study *"Pervasive diversifying selection on the
respiratory syncytial virus G glycoprotein is concentrated in mucin-like flanks, with
episodic sensitivity signals outside CX3C."*

**Summary of findings.** Using a uniform align -> maximum-likelihood phylogeny -> FEL
site-wise selection pipeline, validated on influenza H3N2 haemagglutinin, we find **no
enrichment** of FEL-positive sites in the RSV G central conserved region (CCR) or CX3C motif
in RSV-A, RSV-B, or pooled data (all *p* = 1.0). A powered, threshold-free region-level test
shows the CCR is under significantly stronger purifying selection than the mucin-like domain
in RSV-A (Kruskal-Wallis *p* = 0.027; CCR vs mucin *p* = 0.0066), with a concordant trend in
RSV-B. All A2-mappable FEL-positive RSV sites fall in the C-terminal mucin-like domain.
MEME was added as a secondary episodic-selection sensitivity analysis: it did not identify
episodic diversifying selection in the CX3C motif, but it did flag one non-CX3C broader-CCR
site in each RSV subtype, plus repeated mucin-domain and insertion/unmappable signals.
Structural mapping onto the CX3C/cystine-noose core (PDB 5WN9) confirms the purifying
signature over the functional motif.

## Repository structure

```
.
├── code/            analysis scripts (Python) + PyMOL render script
├── data/            input sequences (FASTA/GenBank) and structure coordinates (PDB 5WN9)
├── results/         alignments, ML trees, FEL/MEME outputs, per-site/region tables (JSON/CSV)
├── figures/         publication figures (PNG)
├── notes/           lab_notebook.txt (full narrative + SHA256 checksums) and accession lists
└── papers/          manuscript draft (Markdown)
```

## Data and code availability

All sequence accessions (`notes/*_accessions.tsv`), alignments (`results/*_aligned*.fasta`),
maximum-likelihood trees (`results/*_iqtree.treefile` and IQ-TREE reports), FEL outputs
(`results/*_fel.json`, `results/*_fel_sites.csv`), MEME outputs (`results/meme/`), analysis scripts (`code/`), figures
(`figures/`), and a complete lab notebook with SHA256 checksums (`notes/lab_notebook.txt`) are
included in this repository. Structure coordinates are from the RCSB PDB, accession **5WN9**
(`data/5WN9.pdb`).

## Requirements

**Python** (see `requirements.txt`):

```
pip install -r requirements.txt
```

- Python ≥ 3.11
- Biopython 1.87, NumPy, SciPy, Matplotlib

**External tools** (invoked outside the Python scripts; versions used in this study):

| Tool | Version | Step |
|---|---|---|
| MAFFT | 7.526 | codon-aware / protein alignment |
| IQ-TREE | 3.1.3 | maximum-likelihood phylogeny + ModelFinder |
| HyPhy (FEL/MEME) | 2.5.100 (Docker `quay.io/biocontainers/hyphy:2.5.100`) | pervasive and episodic site-wise selection |
| PyMOL (open-source) | 3.1.0 | structure rendering |

Best-fit substitution models selected by ModelFinder: **TVM+F+I+R2** (influenza H3 HA),
**TN+F+I+R2** (RSV-A and RSV-B G).

## Reproducing the analysis

The pipeline runs in the following order. Alignment (MAFFT), tree inference (IQ-TREE), and
selection (HyPhy FEL/MEME) are external steps; their inputs and outputs are in `data/` and `results/`.

1. `code/download_sequences.py` — retrieve coding sequences (RSV-A, RSV-B, influenza H3 HA control).
2. `code/extract_cds.py` — extract the G / HA coding regions.
3. `code/codon_align.py` — build codon alignments from protein alignments.
4. **MAFFT** — protein/codon alignment (external).
5. **IQ-TREE** — maximum-likelihood tree with ModelFinder (external).
6. **HyPhy FEL** — per-codon pervasive dN/dS (external).
7. `code/parse_fel.py` — parse FEL JSON and list positively selected sites.
8. `code/map_numbering.py`, `code/map_numbering_rsv.py` — map alignment columns to reference
   coordinates (H3 HA1; RSV-A2 M11486) via the reference added with MAFFT `--add --keeplength`.
9. `code/functional_regions.py` — pre-registered region definitions (CCR 157–198, CX3C 182–186,
   mucin 199–298; H3 antigenic sites A–E). Locked 2026-07-04.
10. `code/enrichment_test.py` — primary permutation enrichment test (10,000 permutations, seed 42).
11. `code/posthoc_mucin_enrichment.py` — exploratory mucin-domain enrichment.
12. `code/region_selection_analysis.py` — powered region-level dN−dS test (Kruskal–Wallis +
    directional Mann–Whitney with rank-biserial effect sizes). Writes
    `results/*_region_selection.json` and per-codon CSVs.
13. **HyPhy MEME** — secondary episodic site-wise sensitivity analysis (external).
14. `code/parse_meme.py` — parse MEME JSON, map RSV sites to RSV-A2 coordinates, and write
    FEL/MEME comparison summaries under `results/meme/`.
15. `code/make_figures.py` — permutation-null histograms, site lollipop, trees, and summary tables.
16. `code/region_selection_figure.py` — per-codon dN−dS distribution figure.
17. `code/structure_figure.py` — Matplotlib structure figure; also writes the B-factor-encoded PDB.
18. `code/render_ccd_pymol.pml` (open-source PyMOL) + `code/compose_structure_figure.py` —
    ray-traced structure panel with colorbar.

Example (powered region-level test):

```
python code/region_selection_analysis.py \
    --fel-json results/rsv_a_fel.json \
    --aln-with-A2 results/rsv_a_prot_aligned_with_A2.fasta \
    --label rsv_a
```

## Key outputs

| Output | File |
|---|---|
| Enrichment summary (primary) | `results/summary_enrichment.md`, `results/summary_enrichment.csv` |
| Region-level selection (powered) | `results/summary_region_selection.csv`, `results/*_region_selection.json` |
| Selected-site list | `results/summary_selected_sites.csv` |
| MEME episodic sensitivity | `results/meme/*_meme.json`, `results/meme/summary_fel_meme_sites.csv`, `results/meme/summary_meme_region_counts.csv` |
| Region-level dN−dS figure | `figures/rsv_region_selection_dnds.png` |
| Structure panel (PyMOL) | `figures/rsv_g_ccd_structure_pymol.png` |
| Full narrative + checksums | `notes/lab_notebook.txt` |
| Manuscript draft | `papers/RSV_G_selection_manuscript.md` |

## License

To be added before publication (e.g., MIT for code, CC-BY for data/figures).

## Citation

If you use this pipeline or data, please cite the manuscript
(`papers/RSV_G_selection_manuscript.md`). Structure coordinates: RCSB PDB 5WN9.
