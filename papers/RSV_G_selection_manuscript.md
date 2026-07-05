# RSV G selection analyses support CX3C/cystine-noose constraint and diversification in mucin-like regions

**Authors:** Zahid [Last name]1*, [Co-authors]

1 [Affiliation]; ORCID: [to add]

\* Correspondence: zahidtzzzz@gmail.com

**Running title:** RSV G selection is partitioned by region

**Keywords:** respiratory syncytial virus; G glycoprotein; central conserved region; CX3C motif; purifying selection; dN/dS; molecular evolution; FEL; MEME

---

> **Reproducibility note.** This manuscript was assembled directly from the analysis pipeline. Numeric results, figures, tables, MEME JSON files, parsing scripts, and SHA256 checksums are traceable to `results/`, `results/meme/`, `figures/`, and `notes/lab_notebook.txt` in the GitHub repository. The region definitions were pre-specified before RSV selected-site positions were inspected. See `papers/reproducibility_note.md` for full details on frozen inputs, HyPhy versions, and checksum verification.

---

## Abstract

**Background:** The RSV attachment (G) glycoprotein is a variable surface antigen and an important target of humoral immunity, while RSV F remains the leading target of licensed and advanced vaccine/immunoprophylaxis approaches. G nevertheless contains a short central conserved region (CCR; RSV-A2 residues 157-198) with a CX3C chemokine-mimic motif (residues 182-186) and a disulfide-bonded cystine noose, making its evolutionary stability a distinct question.

**Methods:** We tested a pre-specified hypothesis that FEL-positive pervasive diversifying selection is enriched in the CCR and/or CX3C motif. A uniform align -> maximum-likelihood phylogeny -> FEL site-wise selection pipeline was validated against influenza H3N2 haemagglutinin. MEME was then run as a secondary episodic site-wise sensitivity analysis using the same codon alignments and trees. Site-wise FEL and MEME p-values are unadjusted and are interpreted as screening-level calls for region summaries.

**Results:** The influenza control behaved as expected: FEL-positive HA1 sites were enriched in antigenic sites A-E (*p* = 0.0081). In RSV G, no mappable FEL-positive sites fell in the CCR or CX3C motif in RSV-A, RSV-B, or the descriptive pooled summary (all one-sided enrichment *p* = 1.0). Because this site-count test has low power with only 3-6 significant sites, this is evidence that the enrichment hypothesis is not supported, not proof that the CCR never experiences diversifying selection. All six A2-mappable FEL-positive RSV sites fell in the C-terminal mucin-like domain, with one additional RSV-A FEL-positive site in the ON1 insertion. A threshold-free companion analysis of per-codon FEL beta-alpha estimates showed stronger CCR purifying constraint in RSV-A (Kruskal-Wallis *p* = 0.0275; CCR vs mucin Mann-Whitney *p* = 0.0066, rank-biserial -0.264) and the same non-significant direction in RSV-B (CCR vs mucin *p* = 0.0742). MEME identified one non-CX3C broader-CCR site in each subtype (RSV-A2 residues 178 and 164), along with mucin-domain and insertion/unmappable sites, but did not identify episodic diversifying selection in the CX3C motif in this dataset.

**Interpretation:** FEL supports concentration of pervasive diversifying selection in mucin-like flanks rather than the CCR/CX3C motif, while MEME adds sparse episodic signals in the broader CCR outside CX3C. Together with structural and antibody studies of the RSV G central conserved domain, these evolutionary data support, but do not by themselves establish, the rationale for monitoring and targeting the CX3C/cystine-noose core.

---

## 1. Introduction

Respiratory syncytial virus (RSV) is a leading cause of lower-respiratory-tract infection in infants, older adults, and immunocompromised individuals. RSV F remains the dominant target of licensed and advanced vaccine and immunoprophylaxis approaches because it is comparatively conserved and carries potent neutralizing epitopes [refs]. RSV G is more variable, but it is also antibody-accessible, immunologically important, and central to attachment biology. That combination makes G a useful test case for whether immune pressure and functional constraint are spatially separated within one glycoprotein.

Two antigenic subgroups, RSV-A and RSV-B, co-circulate and are distinguished largely by their attachment (G) glycoprotein. G is a heavily O-glycosylated type-II membrane protein whose ectodomain is dominated by mucin-like, hypervariable regions flanking a short central conserved region. We use **CCR** for RSV-A2 G residues 157-198, **CX3C motif** for residues 182-186 nested within the CCR, and **CX3C/cystine-noose core** for the ordered structural peptide segment visualized in PDB 5WN9. The CX3C motif mimics the chemokine fractalkine (CX3CL1), engages CX3CR1, and is presented by a cystine-noose scaffold formed by nested disulfide bonds [refs].

The CCR is the epitope for protective monoclonal antibodies and is used in several G-focused vaccine concepts, but those facts do not imply that it is a broad hotspot of adaptive change. If antibody pressure drives adaptive change in G, one possible pattern is enrichment of episodic or pervasive diversifying selection in the CCR, analogous to the antigenic sites of influenza haemagglutinin. A second pattern is stronger functional and structural constraint within the CX3C/cystine-noose core, with diversification accommodated more readily in the mucin-like flanks.

Here we distinguish those patterns with a reproducible molecular-evolution pipeline. The CCR and CX3C region definitions were specified and committed before RSV selected-site positions were inspected. We first validate the FEL pipeline on influenza H3N2 haemagglutinin, where positive selection is expected in antigenic sites, and then apply the same workflow to RSV-A and RSV-B G. Because the pre-specified FEL site-count test is weak when few sites reach significance, we add a threshold-free region-level analysis of per-codon selection summaries, map the result to the ordered CX3C/cystine-noose structure, and use MEME as a secondary episodic-selection sensitivity analysis.

---

## 2. Methods

### 2.1 Sequence data

Frozen coding-sequence inputs for the RSV G gene were used for RSV-A (n = 108 sequences) and RSV-B (n = 115 sequences); influenza A/H3N2 haemagglutinin (HA) coding sequences were used as a positive-control dataset. Accession lists are provided in `notes/rsv_a_nt_accessions.tsv`, `notes/rsv_b_nt_accessions.tsv`, and `notes/flu_h3_ha_nt_accessions.tsv`. RSV FASTA files are archived in `data/` and treated as fixed inputs because the exact original RSV Entrez query text and random seed were not captured at download time. The lab notebook records this limitation explicitly and provides checksums for the frozen FASTA/CDS files. Influenza sampling was scripted in `code/download_sequences.py` with date bins, length filtering, sequence-level deduplication, and seed 42. All RSV G coordinates are reported in **RSV-A2 G protein residue numbering** using GenBank M11486 as the 298-aa reference; H3 HA sites are reported in mature-HA1 (H3) numbering.

### 2.2 Alignment and phylogenetics

Nucleotide coding sequences were codon-aware aligned with **MAFFT v7.526** (translation-guided, back-translated to a codon alignment). Maximum-likelihood phylogenies were inferred with **IQ-TREE v3.1.3** with automated model selection (ModelFinder). The best-fit substitution models were **TVM+F+I+R2** for H3N2 HA and **TN+F+I+R2** for both RSV-A and RSV-B G. Trees were used unrooted for the selection analysis (Supplementary Figures S1-S3). No recombination screen or tree-uncertainty sensitivity analysis was performed.

### 2.3 Site-wise selection

Per-codon pervasive selection was estimated with **FEL (Fixed Effects Likelihood)** in **HyPhy v2.5.100** (run in Docker; image `quay.io/biocontainers/hyphy:2.5.100`). FEL estimates a synonymous rate (alpha = dS) and a nonsynonymous rate (beta = dN) at each codon. A site was called **FEL-positive (pervasively diversifying)** when beta > alpha at unadjusted *p* <= 0.05, and **purifying** when beta < alpha at unadjusted *p* <= 0.05. The per-codon selection metric used throughout is **dN - dS (beta - alpha)**. These p-values are unadjusted codon-wise screening thresholds rather than FDR-controlled discovery claims; individual sites are therefore interpreted cautiously, with emphasis on region-level patterns.

MEME was run as a secondary sensitivity analysis to test for episodic site-wise diversifying selection using the same codon alignments and phylogenies as FEL. MEME used the Universal genetic code, all branches, HyPhy v2.5.100 in the same Docker image, and a reporting threshold of unadjusted *p* <= 0.05. Raw JSON outputs and parsed summaries are stored under `results/meme/`.

### 2.4 Pre-specified functional regions

In RSV-A2 G protein numbering: **N-terminal region** 1-156; **CCR** 157-198; **CX3C motif** 182-186 (nested within the CCR); **mucin-like C-terminal domain** 199-298. For the RSV-A ON1 genotype, the 24-amino-acid C-terminal duplication has no A2-equivalent residue and was labelled "insertion." H3 HA antigenic sites A-E were taken from the canonical Wiley/Wilson definitions [refs]. In RSV-B, A2 mapping covered all CCR/CX3C residues but omitted four distal C-terminal A2 residues; analyses involving the mucin tail should therefore be interpreted with this minor mapping loss in mind.

### 2.5 Primary test: enrichment of selected sites

FEL-significant site positions were mapped to reference coordinates via the reference sequence added to each protein alignment with MAFFT `--add --keeplength`. Enrichment of positively selected sites within a target region was assessed by a **one-sided enrichment permutation test** (10,000 permutations, seed 42): the observed fraction of selected sites falling in the region was compared to a null distribution generated by placing the same number of sites at random over the mapped protein length. Tests were run for the CCR and the CX3C motif in RSV-A and RSV-B. The pooled RSV-A+B rows combine A2-mappable selected sites from subtype-specific analyses as a descriptive summary, not as a joint evolutionary model. Because the permutation test is one-sided for enrichment, *p* = 1.0 supports "not enriched" rather than formal depletion or absence.

### 2.6 Threshold-free companion analysis: region-level selection distribution

Because the primary selected-site test has low power when few sites are significant, we compared the **full per-codon dN - dS distribution** across the pre-specified regions using all mapped codons (script `code/region_selection_analysis.py`). This threshold-free companion analysis uses approximately 300 codons per subtype. An omnibus **Kruskal-Wallis** test across N-terminal / CCR / mucin regions was followed by directional **Mann-Whitney U** tests (a priori: CCR more purifying than mucin; CCR more purifying than the rest of the protein; CX3C more purifying than mucin), each reported with a **rank-biserial** effect size (|r|: 0.1 small, 0.3 medium, 0.5 large). A secondary Fisher exact test compared the fraction of significantly purifying codons between CCR and mucin. This analysis treats FEL per-codon estimates as comparable site-level summaries and does not propagate uncertainty in individual rate estimates or shared phylogenetic/model-fitting uncertainty.

### 2.7 Structural mapping

The experimental structure of RSV-A2 G bound to scFv 2D10 (**PDB 5WN9**) was retrieved from the RCSB. The ordered chain-A coordinates used for coloring span A2 residues 169-189; the deposited peptide includes additional flanking residues that are not represented as the same ordered coordinate span. The modeled sequence is `NFVPCSICSNNPTCWAICKRI`, and disulfides 173-186 and 176-182 were confirmed from coordinates at S-S approximately 2.0 A. Per-codon dN - dS values from the RSV-A FEL analysis were written into the B-factor column and rendered with **open-source PyMOL v3.1.0** (`code/render_ccd_pymol.pml`), with a Matplotlib-composited colorbar and annotation (`code/compose_structure_figure.py`).

### 2.8 Software and reproducibility

Python 3.11 with Biopython 1.87, NumPy, SciPy, and Matplotlib. All inputs, outputs, scripts, and a lab notebook with SHA256 checksums are in the project repository. Every figure is regenerated from `results/*.json` and `results/*.csv` only; MEME summaries are regenerated from raw `results/meme/*_meme.json` by `code/parse_meme.py`; no manuscript values are hand-entered.

---

## 3. Results

### 3.1 The pipeline recovers known positive selection in influenza HA

Applied to influenza H3N2 HA, the FEL pipeline detected five HA1-mappable positively selected codons, **all five (100%) of which fall in antigenic sites A-E** (observed 100% vs 39.4% expected by chance; enrichment 2.5x; permutation **_p_ = 0.0081**; Figure 1, Table 1). The individual sites map to HA1 antigenic sites C (residue 50), A (135), B (157), B (193), and E (261) (Table S1). MEME detected 11 episodic sites in the same influenza dataset, six of which mapped to HA1 antigenic sites. This validates that the align -> tree -> selection workflow can recover immune-driven diversifying selection where it is expected.

### 3.2 FEL-positive selection is not enriched in the RSV G CCR or CX3C motif

In contrast to influenza, **no** FEL-positive site fell within the RSV G CCR or CX3C motif in either subtype or in the descriptive pooled summary (Table 1, Figure 2):

- RSV-A CCR: 0/3 mappable selected sites (observed 0% vs 13.9% expected; *p* = 1.0); CX3C 0/3 (*p* = 1.0).
- RSV-B CCR: 0/3 (*p* = 1.0); CX3C 0/3 (*p* = 1.0).
- RSV-A+B descriptive pooled summary: 0/6 in CCR (0% vs 14.1% expected; *p* = 1.0); 0/6 in CX3C (*p* = 1.0).

The pre-specified FEL enrichment hypothesis was therefore **not supported**: no mappable FEL-positive RSV G sites fell in the CCR or CX3C motif. Because only 3-6 sites were available for the site-count test and the test was one-sided for enrichment, this result should be interpreted as absence of enrichment evidence, not proof that the CCR is incapable of episodic or pervasive diversifying selection.

### 3.3 FEL-positive RSV sites localize to the mucin-like C-terminal domain or ON1 insertion

Every A2-mappable FEL-positive RSV site - six across both subtypes (RSV-A2 positions 237, 247, 274 from RSV-A; 218, 266, 284 from RSV-B) - lies in the **C-terminal mucin-like hypervariable domain** (residues 199-298), i.e., 6/6 versus approximately 34% expected (Figure 2, Table S1). An additional strongly selected RSV-A site (beta = 8.9, *p* = 0.009) falls within the ON1-genotype C-terminal duplication itself, indicating that the duplicated segment is also a locus of diversification. Because the mucin-domain boundary was examined after observing site positions, this localization is reported as **exploratory / post hoc**.

### 3.4 The CCR has lower dN - dS than the mucin domain in RSV-A

The threshold-free region-level analysis, which uses all mapped codons rather than only sites reaching a significance threshold, provides a higher-information companion to the sparse site-count test (Figure 3, Table 2):

- **RSV-A:** the three regions differ in dN - dS (Kruskal-Wallis H = 7.19, *p* = 0.0275). The **CCR is significantly more purifying than the mucin domain** (Mann-Whitney *p* = 0.0066, rank-biserial = -0.264, a small-to-moderate shift) and more purifying than the rest of the protein (*p* = 0.0094, rbc = -0.224).
- **RSV-B:** the same directional pattern holds but does not reach significance (omnibus *p* = 0.252; CCR vs mucin *p* = 0.0742, rbc = -0.154; CCR vs rest *p* = 0.0561).

Thus, region-level dN - dS distributions support stronger CCR purifying constraint in RSV-A and a concordant but non-significant trend in RSV-B.

### 3.5 The purifying signature maps onto the ordered CX3C/cystine-noose structure

Coloring the ordered 5WN9 chain-A coordinates (A2 residues 169-189) by per-codon dN - dS shows the CX3C motif and cystine noose sitting in a predominantly purifying (dN < dS) core, with only two mildly non-negative residues (Figure 4). The two nested disulfides (173-186, 176-182) that lock the fractalkine-mimic motif fall entirely within this constrained segment, providing a structural rationale for the constraint quantified in Section 3.4.

### 3.6 MEME episodic-selection sensitivity analysis

MEME was added after the FEL and region-level analyses as a secondary test for episodic site-wise diversifying selection (Table S2). In RSV-A, MEME identified eight sites at *p* <= 0.05: one A2-mappable CCR site outside the CX3C motif (A2 residue 178; *p* = 0.0438), four A2-mappable mucin-domain sites, one A2-mappable N-terminal site, and two ON1/unmappable insertion sites. In RSV-B, MEME identified seven sites: one A2-mappable CCR site outside the CX3C motif (A2 residue 164; *p* = 0.00053), five mucin-domain sites, and one N-terminal site. MEME did not identify episodic diversifying selection in the CX3C motif in this dataset. The two broader-CCR MEME hits were not FEL-positive and were estimated to involve a small episodic component (approximately one branch under selection in each subtype), so they are interpreted as sensitivity signals rather than a replacement for the FEL/region-level result.

---

## 4. Discussion

The main result is model-specific. FEL-positive pervasive diversification was not enriched in the CCR or CX3C motif and all A2-mappable FEL-positive RSV sites fell in the mucin-like C-terminal domain. MEME, which tests episodic diversifying selection, changes the interpretation for the broader CCR but not for the CX3C motif: episodic-selection calls include one CCR site in each subtype, whereas the CX3C motif has no MEME-positive calls in either subtype. These results are complementary rather than interchangeable.

This pattern is mechanistically coherent. The CX3C motif mimics fractalkine and engages CX3CR1, a function that plausibly imposes purifying selection on the motif and the disulfide scaffold that presents it. By contrast, substitutions and indels may be more readily tolerated in the mucin-like domains, which are long, O-glycosylated, structurally flexible, and known to accommodate insertion/duplication events such as ON1. The influenza HA positive control demonstrates that the pipeline can detect immune-driven diversifying selection when the immunodominant region is not constrained in the same way, while the MEME sensitivity analysis cautions against describing the broader CCR as completely spared from episodic signals.

The subtype asymmetry should be kept explicit. The region-level result is statistically clear only for RSV-A, while RSV-B shows the same direction without reaching significance. This may reflect genuinely weaker constraint, fewer informative substitutions in the RSV-B CCR, genotype or sampling differences, or the minor C-terminal mapping loss disclosed in Methods. The analysis does not justify treating the RSV-B region-level signal as equally supported.

**Implications.** Together with structural and antibody studies of the RSV G central conserved domain, these evolutionary results support the rationale for monitoring and targeting the CX3C/cystine-noose core. They do not by themselves establish vaccine or therapeutic efficacy, and adjacent non-CX3C CCR/noose residues should be described more cautiously because MEME flags non-CX3C episodic sensitivity hits. Surveillance for antigenic change in G should still include the mucin-like domains, where both FEL and MEME identify repeated adaptive signal.

---

## 5. Limitations

1. **Sequence sampling reproducibility.** The exact original RSV download query and random seed were not recorded. RSV FASTA/accession files are therefore frozen inputs with checksums, not fully reproducible downloads. A submission-ready version should add database date, exact query, inclusion/exclusion rules, genotype/year/geography summaries, deduplication rules, ambiguous-base handling, and accession-list checksums in a supplement.
2. **Power of the primary selected-site test.** With only 3-6 FEL-significant sites per subtype, the site-count enrichment test has limited power. A one-sided enrichment *p* = 1.0 should be read as "not enriched," not as proof of absence or formal depletion.
3. **Unadjusted site-wise thresholds.** FEL and MEME site-wise p-values are unadjusted. They are useful as screening thresholds for region summaries, but individual selected sites may be sensitive to multiple testing, alignment, and model assumptions.
4. **Threshold-free companion analysis assumptions.** The region-level analysis uses all codons, but it treats FEL point estimates as comparable site-level summaries and does not propagate uncertainty in individual rate estimates, tree uncertainty, or shared phylogenetic/model-fitting uncertainty.
5. **Recombination and phylogenetic sensitivity.** No recombination screen, temporal/genotype clustering analysis, duplicate/near-duplicate sensitivity analysis, or alternative-tree sensitivity analysis was performed.
6. **Subtype dependence.** The purifying-constraint result is statistically clear only for RSV-A; RSV-B shows the same direction as a non-significant trend.
7. **Cross-subtype coordinate mapping.** RSV-B G was mapped onto RSV-A2 numbering. Adding the A2 reference via MAFFT `--keeplength` dropped four distal C-terminal A2 residues (a2_ref_len 294 vs 298), affecting a small part of the mucin tail but not the CCR/CX3C boundaries.
8. **Structural coverage.** No full-length or mucin-domain G structure is included here; the ordered 5WN9 coordinates used in Figure 4 span A2 residues 169-189. The structural figure therefore illustrates constraint on the CX3C/cystine-noose core, not the entire CCR or flanks.
9. **Post hoc mucin localization.** The concentration of selected sites in the mucin domain used a boundary emphasized after observing site positions and is reported as exploratory.
10. **MEME sensitivity.** MEME detects episodic site-wise selection and found one non-CX3C CCR site in each RSV subtype. These calls are not FEL-positive, are branch-sparse, and may be sensitive to alignment, duplicate sequences, zero-length branch handling, and A2 mapping; they should be treated as sensitivity signals rather than proof that the whole CCR is adaptively labile.

---

## 6. Conclusion

FEL-detected pervasive diversifying selection on the RSV G glycoprotein is largely outside the CCR/CX3C motif and is concentrated among mucin-like or insertion/unmappable positions in this dataset. The threshold-free region-level analysis supports stronger CCR purifying constraint in RSV-A and a concordant but non-significant trend in RSV-B. MEME did not identify episodic diversifying selection in the CX3C motif in this dataset, but it did flag one non-CX3C broader-CCR site in each subtype. These results support the CX3C/cystine-noose core as a comparatively constrained region to monitor and study, while keeping vaccine and therapeutic claims dependent on independent experimental evidence.

---

## Figures

**Figure 1. Pipeline validation on influenza H3N2 HA.** Permutation-null distribution of the fraction of positively selected HA sites expected in antigenic sites A-E; the observed fraction (red line, 100%) is far in the upper tail (*p* = 0.0081). Source figure: `figures/flu_null_hist.png`.

**Figure 2. FEL-positive RSV G sites relative to the CCR and CX3C motif.** (a) One-sided enrichment permutation null for the CCR in the descriptive pooled RSV summary; observed enrichment = 0 (*p* = 1.0). (b) Lollipop plot of every FEL-positive site along the RSV-A2 G coordinate axis, with the CCR and CX3C motif shaded; all A2-mappable sites fall in the C-terminal mucin-like domain. Source figures: `figures/rsv_combined_ccr_null_hist.png`, `figures/rsv_selected_sites_lollipop.png`.

**Figure 3. Region-level selection pressure (threshold-free companion analysis).** Per-codon dN - dS distributions by region (N-terminal / CCR / mucin) for RSV-A and RSV-B, with Kruskal-Wallis and CCR-vs-mucin Mann-Whitney statistics annotated. The CCR box is shifted toward purifying selection (below zero), significantly so in RSV-A. Source figure: `figures/rsv_region_selection_dnds.png`.

**Figure 4. The ordered CX3C/cystine-noose core is purifying.** Ray-traced structure of the ordered RSV-A2 G CCD coordinates in PDB 5WN9 (A2 residues 169-189) with the backbone colored by per-codon dN - dS (blue = purifying, red = diversifying); the two nested disulfides (173-186, 176-182) are shown in gold and the CX3C motif is labelled. Source figure: `figures/rsv_g_ccd_structure_pymol.png`.

**Supplementary figures.** Maximum-likelihood phylogenies (`figures/flu_h3_ha_tree.png`, `figures/rsv_a_tree.png`, `figures/rsv_b_tree.png`); per-subtype CCR and CX3C permutation-null histograms (`figures/rsv_a_ccr_null_hist.png`, `figures/rsv_a_cx3c_null_hist.png`, `figures/rsv_b_ccr_null_hist.png`, `figures/rsv_b_cx3c_null_hist.png`, `figures/rsv_combined_cx3c_null_hist.png`).

---

## Tables

**Table 1. Enrichment of FEL-positive sites (one-sided enrichment permutation test, 10,000 permutations, seed 42).**

| Target | Region | n | hits | obs % | exp % | enrichment | p-value | sig |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| flu | H3 antigenic sites A-E | 5 | 5 | 100.0 | 39.43 | 2.54x | 0.0081 | **YES** |
| RSV-A | CCR 157-198 | 3 | 0 | 0.0 | 13.9 | 0x | 1.0 | no |
| RSV-A | CX3C 182-186 | 3 | 0 | 0.0 | 1.58 | 0x | 1.0 | no |
| RSV-B | CCR 157-198 | 3 | 0 | 0.0 | 13.9 | 0x | 1.0 | no |
| RSV-B | CX3C 182-186 | 3 | 0 | 0.0 | 1.58 | 0x | 1.0 | no |
| RSV-A+B descriptive pooled summary | CCR 157-198 | 6 | 0 | 0.0 | 14.1 | 0x | 1.0 | no |
| RSV-A+B descriptive pooled summary | CX3C 182-186 | 6 | 0 | 0.0 | 1.66 | 0x | 1.0 | no |

Source data: `results/summary_enrichment.md` and `results/summary_enrichment.csv`.

**Table 2. Region-level selection (per-codon dN - dS; threshold-free companion analysis).**

| Target | Comparison | Test | Statistic | Effect size | p-value | sig |
|---|---|---|---|---|---:|:---:|
| RSV-A | N-term/CCR/mucin (omnibus) | Kruskal-Wallis | H = 7.19 | epsilon2 = 0.018 | 0.0275 | **YES** |
| RSV-A | CCR < mucin (more purifying) | Mann-Whitney U | U = 1546 | rbc = -0.264 | 0.0066 | **YES** |
| RSV-A | CCR < rest of protein | Mann-Whitney U | U = 4172 | rbc = -0.224 | 0.0094 | **YES** |
| RSV-B | N-term/CCR/mucin (omnibus) | Kruskal-Wallis | H = 2.76 | epsilon2 = 0.003 | 0.252 | no |
| RSV-B | CCR < mucin (more purifying) | Mann-Whitney U | U = 1705 | rbc = -0.154 | 0.0742 | no |
| RSV-B | CCR < rest of protein | Mann-Whitney U | U = 4494 | rbc = -0.151 | 0.0561 | no |

Source data: `results/summary_region_selection.csv`. Negative rank-biserial (rbc) indicates lower CCR dN - dS than the comparator region.

**Table S1.** Full list of FEL-positive sites (flu + RSV-A + RSV-B) with reference coordinates and beta estimates: `results/summary_selected_sites.csv`.

**Table S2.** MEME episodic-selection sensitivity summary and FEL/MEME positive-site comparison: `results/meme/summary_meme_region_counts.csv` and `results/meme/summary_fel_meme_sites.csv`.

---

## Data and code availability

All sequence accessions, frozen FASTA inputs, alignments, trees, FEL outputs, MEME outputs, analysis scripts, figures, and the lab notebook with SHA256 checksums are provided in the project repository: <https://github.com/zahid-bio/rsv-g-selection-analysis> (Zenodo release: [to add DOI after preprint posting]). Structure coordinates are from the RCSB PDB (accession 5WN9). Code is licensed under the MIT License; all manuscript text, figures, and data are licensed under CC-BY 4.0.

## Author contributions

[To complete.]

## Funding

[To complete.]

## Acknowledgements

None. All analyses were performed independently using published algorithms and publicly available data.

## Conflicts of interest

The authors declare no competing interests.

---

## References

> **These references are a canonical scaffold and must be checked against the primary literature (authors, year, volume, pages, DOI) before submission.**

1. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7. *Mol Biol Evol.* 2013;30(4):772-780.
2. Minh BQ, et al. IQ-TREE 2: New models and efficient methods for phylogenetic inference. *Mol Biol Evol.* 2020;37(5):1530-1534.
3. Kalyaanamoorthy S, et al. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nat Methods.* 2017;14:587-589.
4. Kosakovsky Pond SL, Frost SDW. Not so different after all: a comparison of methods for detecting amino acid sites under selection. *Mol Biol Evol.* 2005;22(5):1208-1222. (FEL)
5. Kosakovsky Pond SL, et al. HyPhy: hypothesis testing using phylogenies. *Bioinformatics.* 2005;21(5):676-679.
6. Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics.* 2009;25(11):1422-1423.
7. Tripp RA, et al. CX3C chemokine mimicry by respiratory syncytial virus G glycoprotein. *Nat Immunol.* 2001;2(8):732-738.
8. [RSV G central conserved domain structure / scFv 2D10; PDB 5WN9] - verify primary citation.
9. Wiley DC, Wilson IA, Skehel JJ. Structural identification of the antibody-binding sites of Hong Kong influenza haemagglutinin and their involvement in antigenic variation. *Nature.* 1981;289:373-378. (H3 antigenic sites)
10. The PyMOL Molecular Graphics System, open-source build, v3.1.0. Schrodinger, LLC.
11. Murrell B, et al. Detecting individual sites subject to episodic diversifying selection. *PLoS Genet.* 2012;8(7):e1002764. (MEME)
12. [RSV F as the leading antibody/vaccine target; Nature Communications 2019] - verify citation details.
13. [RSV G central conserved domain vaccine context; npj Vaccines 2022] - verify citation details.
14. [RSV G functional-domain review, including CCR/CX3C definitions; Viruses 2021] - verify citation details.
15. Fedechkin SO, et al. [RSV G central conserved domain immunogen/antibody study]. *Sci Immunol.* 2018 - verify citation details.
16. [Protective antigenic sites outside the RSV G central conserved domain; PLoS Pathog.] - verify citation details.

---

*Manuscript draft revised 2026-07-04 after critical review of `main (5).pdf`. Verify all references and fill placeholder author/affiliation/funding fields before submission.*
