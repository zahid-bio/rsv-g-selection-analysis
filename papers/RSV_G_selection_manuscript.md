# Pervasive diversifying selection on the respiratory syncytial virus G glycoprotein is concentrated in mucin-like flanks, with episodic sensitivity signals outside CX3C

**Authors:** Zahid [Last name]¹*, [Co-authors]

¹ [Affiliation]

\* Correspondence: zahidtzzzz@gmail.com

**Running title:** RSV G selection is partitioned by region

**Keywords:** respiratory syncytial virus; G glycoprotein; central conserved region; CX3C motif; purifying selection; dN/dS; molecular evolution; FEL; MEME

---

> **DRAFT — status note.** This manuscript was assembled directly from the finished analysis pipeline (`U6_project/`). All numeric results, figures, and tables are traceable to `results/`, `figures/`, and `notes/lab_notebook.txt` (with SHA256 hashes). **The reference list is a scaffold of canonical sources and MUST be verified (authors, years, volumes, DOIs) before submission.** Author names, affiliations, and funding are placeholders.

---

## Abstract

The attachment (G) glycoprotein of respiratory syncytial virus (RSV) is the principal target of the antibody response and the most variable protein in the RSV genome, yet it carries a short central conserved region (CCR, residues 157–198 in the RSV-A2 reference) that harbours a CX3C chemokine-mimic motif and a disulfide-bonded "cystine noose." We tested a pre-registered hypothesis that diversifying selection is enriched in the CCR and/or the CX3C motif. Using a uniform align → maximum-likelihood phylogeny → FEL site-wise selection pipeline, validated against influenza H3N2 haemagglutinin (FEL-positive sites significantly enriched in antigenic sites A–E; *p* = 0.0081), we found **no enrichment** of FEL-positive sites in the RSV G CCR or CX3C motif in RSV-A, RSV-B, or the pooled data (all *p* = 1.0). Because that site-count test is underpowered (3–6 significant sites per subtype), we added a powered, threshold-free companion analysis comparing the per-codon dN−dS distribution across regions. In RSV-A the CCR is under significantly stronger purifying selection than the C-terminal mucin-like domain (Kruskal–Wallis *p* = 0.027; CCR vs mucin Mann–Whitney *p* = 0.0066, rank-biserial −0.26); in RSV-B the same direction is present but does not reach significance (*p* = 0.074). All FEL-positive, mappable RSV sites (6/6 across both subtypes) fall in the C-terminal mucin-like hypervariable domain. MEME, added as an episodic site-wise sensitivity analysis, identified one broader-CCR site in each RSV subtype (RSV-A2 residues 178 and 164), neither in the CX3C motif, alongside mucin-domain and insertion/unmappable sites. Mapping per-codon dN−dS onto the experimental structure of the CX3C/cystine-noose core (PDB 5WN9) confirms a purifying-selection signature over the CX3C motif under FEL. We conclude that pervasive FEL-detected diversification of RSV G is concentrated in the mucin-like flanks, while MEME adds limited branch-sparse episodic CCR signals outside CX3C; the CX3C motif remains comparatively constrained, but the broader CCR should not be described as uniformly free of episodic signal.

---

## 1. Introduction

Respiratory syncytial virus (RSV) is a leading cause of lower-respiratory-tract infection in infants, older adults, and immunocompromised individuals. Two antigenic subgroups, RSV-A and RSV-B, co-circulate and are distinguished largely by their attachment (G) glycoprotein. G is a heavily O-glycosylated type-II membrane protein whose ectodomain is dominated by two mucin-like, hypervariable regions flanking a short **central conserved region (CCR)**. The CCR (residues ~157–198 in the RSV-A2 reference, M11486) is the least variable segment of an otherwise highly variable protein and contains two functionally critical features: a **CX3C chemokine motif** (residues 182–186, Cys-X-X-X-Cys) that mimics the chemokine fractalkine (CX3CL1) and engages its receptor CX3CR1, and a **cystine noose** formed by two nested disulfide bonds that structurally locks the motif [refs].

The G protein is the dominant target of the neutralizing and non-neutralizing antibody response, and the CCR in particular is the epitope for a class of protective monoclonal antibodies and a component of several G-based vaccine candidates [refs]. This creates a clear, testable evolutionary prediction. If antibody pressure drives adaptive change in G, one of two patterns should hold: either (i) the immunodominant CCR is a hotspot of episodic positive (diversifying) selection, as is seen in the antigenic sites of influenza haemagglutinin; or (ii) functional and structural constraint on the CX3C motif and cystine noose overrides immune pressure, confining diversifying selection to the mucin-like flanks and leaving the core under purifying selection.

Here we discriminate between these hypotheses with a uniform, reproducible molecular-evolution pipeline. We **pre-registered** the primary hypothesis — enrichment of positively selected sites in the CCR and CX3C motif — and its region definitions (locked 2026-07-04, before the enrichment tests were run). We first validate the FEL pipeline on influenza H3N2 haemagglutinin, where the expected answer (positive selection concentrated in antigenic sites) is known, and then apply it to RSV-A and RSV-B G. Recognising that the pre-registered site-count test is statistically weak when few sites reach significance, we complement it with a powered, threshold-free analysis of the selection-pressure distribution across regions, anchor the result to the experimental structure of the CX3C core, and add MEME as a secondary episodic-selection sensitivity analysis.

---

## 2. Methods

### 2.1 Sequence data
Coding sequences of the RSV G gene were obtained for RSV-A (n = 108 sequences) and RSV-B (n = 115 sequences); influenza A/H3N2 haemagglutinin (HA) coding sequences were obtained as a positive-control dataset. Accession lists are provided in `notes/rsv_a_nt_accessions.tsv`, `notes/rsv_b_nt_accessions.tsv`, and `notes/flu_h3_ha_nt_accessions.tsv`. All coordinates for RSV G are reported in **RSV-A2 (GenBank M11486) mature-G numbering** (298 residues); H3 HA sites are reported in mature-HA1 (H3) numbering.

### 2.2 Alignment and phylogenetics
Nucleotide coding sequences were codon-aware aligned with **MAFFT v7.526** (translation-guided, back-translated to a codon alignment). Maximum-likelihood phylogenies were inferred with **IQ-TREE v3.1.3** with automated model selection (ModelFinder). The best-fit substitution models were **TVM+F+I+R2** for H3N2 HA and **TN+F+I+R2** for both RSV-A and RSV-B G. Trees were used unrooted for the selection analysis (Figures S1–S3).

### 2.3 Site-wise selection
Per-codon pervasive selection was estimated with **FEL (Fixed Effects Likelihood)** in **HyPhy v2.5.100** (run in Docker; image `quay.io/biocontainers/hyphy:2.5.100`). FEL estimates a synonymous rate (α = dS) and a non-synonymous rate (β = dN) at each codon. A site was called **FEL-positive (pervasively diversifying)** when β > α at *p* ≤ 0.05, and **purifying** when β < α at *p* ≤ 0.05. The per-codon selection metric used throughout is **dN − dS (β − α)**.

MEME was run as a secondary sensitivity analysis to test for episodic site-wise diversifying selection using the same codon alignments and phylogenies as FEL. MEME used the Universal genetic code, all branches, HyPhy v2.5.100 in the same Docker image, and a reporting threshold of *p* ≤ 0.05; raw JSON outputs and parsed summaries are stored under `results/meme/`.

### 2.4 Pre-registered functional regions (locked 2026-07-04)
In RSV-A2 mature-G numbering: **N-terminal region** 1–156; **CCR** 157–198; **CX3C motif** 182–186 (nested within the CCR); **mucin-like C-terminal domain** 199–298. For the RSV-A ON1 genotype, the 24-amino-acid C-terminal duplication has no A2-equivalent residue and was labelled "insertion." H3 HA antigenic sites A–E were taken from the canonical Wiley/Wilson definitions [refs].

### 2.5 Primary test: enrichment of selected sites (pre-registered)
FEL-significant site positions were mapped to reference coordinates via the reference sequence added to each protein alignment with MAFFT `--add --keeplength`. Enrichment of positively selected sites within a target region was assessed by a **one-sided permutation test** (10,000 permutations, seed 42): the observed fraction of selected sites falling in the region was compared to a null distribution generated by placing the same number of sites at random over the mapped protein length. Tests were run for the CCR and the CX3C motif in RSV-A, RSV-B, and the pooled ("combined") data.

### 2.6 Powered companion test: region-level selection distribution
Because the primary test has low power when few sites are significant, we compared the **full per-codon dN−dS distribution** across the pre-registered regions using all mapped codons (script `code/region_selection_analysis.py`). This is threshold-free and uses ~300 codons per subtype. An omnibus **Kruskal–Wallis** test across N-terminal / CCR / mucin regions was followed by directional **Mann–Whitney U** tests (a priori: CCR more purifying than mucin; CCR more purifying than the rest of the protein; CX3C more purifying than mucin), each reported with a **rank-biserial** effect size (|r|: 0.1 small, 0.3 medium, 0.5 large). A secondary Fisher exact test compared the fraction of significantly purifying codons between CCR and mucin. Nonparametric tests were chosen a priori given the heavy-tailed, non-normal dN−dS distributions.

### 2.7 Structural mapping
The experimental structure of the RSV-A2 G CCD peptide (**PDB 5WN9**, chain A, residues 169–189, bound to scFv 2D10) was retrieved from the RCSB (identity and residue range verified programmatically: chain A source organism *Human respiratory syncytial virus A2*, sequence `NFVPCSICSNNPTCWAICKRI`; disulfides 173–186 and 176–182 confirmed from coordinates at S–S ≈ 2.0 Å). Per-codon dN−dS values from the RSV-A FEL analysis were written into the B-factor column and rendered with **open-source PyMOL v3.1.0** (`code/render_ccd_pymol.pml`), with a Matplotlib-composited colorbar and annotation (`code/compose_structure_figure.py`).

### 2.8 Software and reproducibility
Python 3.11 with Biopython 1.87, NumPy, SciPy, and Matplotlib. All inputs, outputs, scripts, and a full lab notebook with SHA256 checksums are in the project repository. Every figure is regenerated from `results/*.json` and `results/*.csv` only; MEME summaries are regenerated from raw `results/meme/*_meme.json` by `code/parse_meme.py`; no values are hand-entered.

---

## 3. Results

### 3.1 The pipeline recovers known positive selection in influenza HA (validation)
Applied to influenza H3N2 HA, the FEL pipeline detected five HA1-mappable positively selected codons, **all five (100%) of which fall in antigenic sites A–E** (observed 100% vs 39.4% expected by chance; enrichment 2.5×; permutation **_p_ = 0.0081**; Figure 1, Table 1). The individual sites map to HA1 antigenic sites C (residue 50), A (135), B (157), B (193), and E (261) (Table S1). MEME detected 11 episodic sites in the same influenza dataset, six of which mapped to HA1 antigenic sites. This confirms that the align → tree → selection workflow detects immune-driven diversifying selection where it is known to exist.

### 3.2 FEL-positive selection is not enriched in the RSV G CCR or CX3C motif (primary result)
In contrast to influenza, **no** FEL-positive site fell within the RSV G CCR or CX3C motif in either subtype or in the pooled data (Table 1, Figure 2):

- RSV-A CCR: 0/3 mappable selected sites (observed 0% vs 13.9% expected; *p* = 1.0); CX3C 0/3 (*p* = 1.0).
- RSV-B CCR: 0/3 (*p* = 1.0); CX3C 0/3 (*p* = 1.0).
- Combined: 0/6 in CCR (0% vs 14.1% expected; *p* = 1.0); 0/6 in CX3C (*p* = 1.0).

The pre-registered FEL enrichment hypothesis is therefore **rejected**: FEL-positive pervasive diversifying selection is not concentrated in the conserved core. (We note this site-count test is underpowered given only 3–6 significant sites per subtype; §3.4 addresses this directly.)

### 3.3 All FEL-positive RSV sites localize to the mucin-like C-terminal domain or ON1 insertion
Every A2-mappable FEL-positive RSV site — six across both subtypes (RSV-A2 positions 237, 247, 274 from RSV-A; 218, 266, 284 from RSV-B) — lies in the **C-terminal mucin-like hypervariable domain** (residues 199–298), i.e., 6/6 versus ~34% expected (Figure 2, Table S1). An additional strongly selected RSV-A site (β = 8.9, *p* = 0.009) falls within the ON1-genotype C-terminal duplication itself, indicating that the duplicated segment is also a locus of diversification. Because the mucin-domain boundary was examined after observing site positions, this localization is reported as **exploratory / post-hoc**.

### 3.4 The CCR is under stronger purifying selection than the mucin domain (powered result)
The threshold-free, region-level analysis — which uses all ~300 codons rather than the handful reaching significance — provides the statistically powered complement to the negative site-count test (Figure 3, Table 2):

- **RSV-A:** the three regions differ in dN−dS (Kruskal–Wallis H = 7.19, *p* = 0.027). The **CCR is significantly more purifying than the mucin domain** (Mann–Whitney *p* = 0.0066, rank-biserial = −0.26, a medium effect) and more purifying than the rest of the protein (*p* = 0.0094, rbc = −0.22).
- **RSV-B:** the same directional pattern holds but does not reach significance (omnibus *p* = 0.25; CCR vs mucin *p* = 0.074, rbc = −0.15; CCR vs rest *p* = 0.056).

Thus the conserved core is not merely devoid of positive selection — in RSV-A it is actively constrained by purifying selection relative to the flanks, with RSV-B showing a concordant but weaker trend.

### 3.5 The purifying signature maps onto the CX3C/cystine-noose structure
Colouring the experimental CCD structure (PDB 5WN9, residues 169–189) by per-codon dN−dS shows the CX3C motif and cystine noose sitting in a predominantly purifying (dN < dS) core, with only two mildly non-negative residues (Figure 4). The two nested disulfides (173–186, 176–182) that lock the fractalkine-mimic motif fall entirely within this constrained segment, providing a structural rationale for the constraint quantified in §3.4.

### 3.6 MEME episodic-selection sensitivity analysis
MEME was added after the FEL and region-level analyses as a secondary test for episodic site-wise diversifying selection (Table S2). In RSV-A, MEME identified eight sites at *p* ≤ 0.05: one A2-mappable CCR site outside the CX3C motif (A2 residue 178; *p* = 0.0438), four A2-mappable mucin-domain sites, one A2-mappable N-terminal site, and two ON1/unmappable insertion sites. In RSV-B, MEME identified seven sites: one A2-mappable CCR site outside the CX3C motif (A2 residue 164; *p* = 0.00053), five mucin-domain sites, and one N-terminal site. MEME did not identify episodic diversifying selection in the CX3C motif in this dataset. The two broader-CCR MEME hits were not FEL-positive and were estimated to involve a small episodic component (approximately one branch under selection in each subtype), so they are interpreted as sensitivity signals rather than a replacement for the FEL/region-level result.

---

## 4. Discussion

Our central finding is that adaptive diversification of the RSV G glycoprotein is **spatially partitioned**, but the partition depends on the selection model. FEL-positive pervasive diversification is excluded from the CCR and CX3C motif and concentrated in the mucin-like flanking domains, while the powered companion analysis shows that the CCR is shifted toward purifying constraint. MEME adds a more sensitive branch-site view: it detects one episodic broader-CCR site in each subtype, neither in CX3C, and most MEME-positive RSV sites still lie in the mucin domain or insertion/unmappable positions. These results are complementary rather than interchangeable.

This pattern is mechanistically coherent. The CX3C motif mimics fractalkine and engages CX3CR1, a function that plausibly imposes strong purifying selection on the motif and the disulfide scaffold that presents it; a diversifying substitution here would risk abolishing receptor engagement and disrupting the noose. Immune escape, by contrast, is cheaply achieved in the mucin-like domains, which are long, O-glycosylated, structurally disordered, and tolerant of substitution, insertion (e.g., the ON1 24-amino-acid duplication), and deletion. The influenza HA positive control demonstrates that the pipeline readily detects immune-driven diversifying selection when the immunodominant region is *not* structurally constrained, while the MEME sensitivity analysis cautions against overstating the broader CCR as entirely free of episodic signals.

The subtype asymmetry (significant constraint in RSV-A, a concordant trend in RSV-B) is worth stating plainly rather than glossing. It may reflect genuinely weaker constraint, smaller effective sample of informative substitutions in the RSV-B CCR, or the minor mapping loss noted below; we do not over-interpret it.

**Implications.** The absence of FEL-positive and MEME-positive sites in the CX3C motif strengthens the rationale for CX3C-motif-directed vaccines and monoclonal antibodies: the functional motif is constrained relative to the variable flanks. Adjacent noose/CCR residues should be described more cautiously because MEME flags non-CX3C episodic sensitivity hits. Surveillance for antigenic change in G should still focus heavily on the mucin-like domains, where both FEL and MEME identify repeated adaptive signal.

---

## 5. Limitations

1. **Power of the primary test.** With only 3–6 FEL-significant sites per subtype, the pre-registered site-count enrichment test has limited power; a null result from it alone would be weak. This is precisely why the region-level distribution test (§3.4), which uses all codons and no significance threshold, is the load-bearing positive analysis.
2. **Subtype dependence.** The purifying-constraint result is statistically clear only for RSV-A; RSV-B shows the same direction as a non-significant trend.
3. **Cross-subtype coordinate mapping.** RSV-B G was mapped onto RSV-A2 numbering (the CCR is conserved across subtypes). Adding the A2 reference via MAFFT `--keeplength` dropped four C-terminal A2 residues in the RSV-B alignment (a2_ref_len 294 vs 298), affecting 4 of ~96 mucin codons; this does not alter the RSV-B conclusion (null regardless) but is disclosed for completeness.
4. **Structural coverage.** No full-length or mucin-domain G structure exists; the crystallized peptide (PDB 5WN9) spans only the CX3C/cystine-noose core (169–189). The structural figure therefore illustrates constraint on the functional core, not the entire CCR or the flanks.
5. **Post-hoc localization.** The concentration of selected sites in the mucin domain (§3.3) used a boundary chosen after observing site positions and is reported as exploratory.
6. **MEME sensitivity.** MEME detects episodic site-wise selection and found one non-CX3C CCR site in each RSV subtype. These calls are not FEL-positive, are branch-sparse, and may be sensitive to alignment, duplicate sequences, zero-length branch handling, and A2 mapping; they should be treated as sensitivity signals rather than proof that the whole CCR is adaptively labile.

---

## 6. Conclusion

FEL-detected pervasive diversifying selection on the RSV G glycoprotein is concentrated in the mucin-like flanks and ON1 insertion, while the CCR is shifted toward purifying constraint (significant in RSV-A, trending in RSV-B). MEME does not identify episodic diversifying selection in the CX3C motif in this dataset, but it does flag one non-CX3C broader-CCR site in each subtype. These results support the CX3C motif as a comparatively constrained vaccine and antibody target, while arguing that antigenic surveillance should remain focused on the variable mucin-like domains and that adjacent CCR/noose positions should be interpreted with model-specific caution.

---

## Figures

**Figure 1. Pipeline validation on influenza H3N2 HA.** Permutation-null distribution of the fraction of positively selected HA sites expected in antigenic sites A–E; the observed fraction (red line, 100%) is far in the upper tail (*p* = 0.0081). `figures/flu_null_hist.png`

**Figure 2. FEL-positive RSV G sites avoid the conserved core.** (a) Permutation-null for the CCR in the pooled RSV data; observed enrichment = 0 (*p* = 1.0). (b) Lollipop plot of every FEL-positive site along the RSV-A2 G coordinate axis, with the CCR and CX3C motif shaded; all A2-mappable sites fall in the C-terminal mucin-like domain. `figures/rsv_combined_ccr_null_hist.png`, `figures/rsv_selected_sites_lollipop.png`

**Figure 3. Region-level selection pressure (powered analysis).** Per-codon dN−dS distributions by region (N-terminal / CCR / mucin) for RSV-A and RSV-B, with Kruskal–Wallis and CCR-vs-mucin Mann–Whitney statistics annotated. The CCR box is shifted toward purifying selection (below zero), significantly so in RSV-A. `figures/rsv_region_selection_dnds.png`

**Figure 4. The CX3C/cystine-noose core is purifying.** Ray-traced structure of the RSV-A2 G CCD (PDB 5WN9, residues 169–189) with the backbone coloured by per-codon dN−dS (blue = purifying, red = diversifying); the two nested disulfides (173–186, 176–182) are shown in gold and the CX3C motif is labelled. `figures/rsv_g_ccd_structure_pymol.png`

**Supplementary figures.** Maximum-likelihood phylogenies (`figures/flu_h3_ha_tree.png`, `figures/rsv_a_tree.png`, `figures/rsv_b_tree.png`); per-subtype CCR and CX3C permutation-null histograms (`figures/rsv_a_ccr_null_hist.png`, `figures/rsv_a_cx3c_null_hist.png`, `figures/rsv_b_ccr_null_hist.png`, `figures/rsv_b_cx3c_null_hist.png`, `figures/rsv_combined_cx3c_null_hist.png`).

---

## Tables

**Table 1. Enrichment of FEL-positive sites (permutation test, 10,000 permutations, seed 42, one-sided).**

| Target | Region | n | hits | obs % | exp % | enrichment | p-value | sig |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| flu | H3 antigenic sites A–E | 5 | 5 | 100.0 | 39.43 | 2.54× | 0.0081 | **YES** |
| RSV-A | CCR 157–198 | 3 | 0 | 0.0 | 13.9 | 0× | 1.0 | no |
| RSV-A | CX3C 182–186 | 3 | 0 | 0.0 | 1.58 | 0× | 1.0 | no |
| RSV-B | CCR 157–198 | 3 | 0 | 0.0 | 13.9 | 0× | 1.0 | no |
| RSV-B | CX3C 182–186 | 3 | 0 | 0.0 | 1.58 | 0× | 1.0 | no |
| RSV-A+B | CCR 157–198 | 6 | 0 | 0.0 | 14.1 | 0× | 1.0 | no |
| RSV-A+B | CX3C 182–186 | 6 | 0 | 0.0 | 1.66 | 0× | 1.0 | no |

Source: `results/summary_enrichment.md`.

**Table 2. Region-level selection (per-codon dN−dS; powered companion).**

| Target | Comparison | Test | Statistic | Effect size | p-value | sig |
|---|---|---|---|---|---:|:---:|
| RSV-A | N-term/CCR/mucin (omnibus) | Kruskal–Wallis | H = 7.19 | ε² = 0.018 | 0.0275 | **YES** |
| RSV-A | CCR < mucin (more purifying) | Mann–Whitney U | U = 1546 | rbc = −0.264 | 0.0066 | **YES** |
| RSV-A | CCR < rest of protein | Mann–Whitney U | U = 4172 | rbc = −0.224 | 0.0094 | **YES** |
| RSV-B | N-term/CCR/mucin (omnibus) | Kruskal–Wallis | H = 2.76 | ε² = 0.003 | 0.252 | no |
| RSV-B | CCR < mucin (more purifying) | Mann–Whitney U | U = 1705 | rbc = −0.154 | 0.0742 | no |
| RSV-B | CCR < rest of protein | Mann–Whitney U | U = 4494 | rbc = −0.151 | 0.0561 | no |

Source: `results/summary_region_selection.csv`. Negative rank-biserial (rbc) = CCR under stronger purifying selection.

**Table S1.** Full list of FEL-positive sites (flu + RSV-A + RSV-B) with reference coordinates and β estimates: `results/summary_selected_sites.csv`.

**Table S2.** MEME episodic-selection sensitivity summary and FEL/MEME positive-site comparison: `results/meme/summary_meme_region_counts.csv` and `results/meme/summary_fel_meme_sites.csv`.

---

## Data and code availability
All sequence accessions, alignments, trees, FEL outputs, MEME outputs, analysis scripts, figures, and a complete lab notebook with SHA256 checksums are provided in the project repository (`U6_project/`). Structure coordinates are from the RCSB PDB (accession 5WN9).

## Author contributions
[To complete.]

## Funding
[To complete.]

## Acknowledgements
[To complete.]

## Conflicts of interest
The authors declare no competing interests.

---

## References

> **These references are a canonical scaffold and must be checked against the primary literature (authors, year, volume, pages, DOI) before submission.**

1. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7. *Mol Biol Evol.* 2013;30(4):772–780.
2. Minh BQ, et al. IQ-TREE 2: New models and efficient methods for phylogenetic inference. *Mol Biol Evol.* 2020;37(5):1530–1534.
3. Kalyaanamoorthy S, et al. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nat Methods.* 2017;14:587–589.
4. Kosakovsky Pond SL, Frost SDW. Not so different after all: a comparison of methods for detecting amino acid sites under selection. *Mol Biol Evol.* 2005;22(5):1208–1222. (FEL)
5. Kosakovsky Pond SL, et al. HyPhy: hypothesis testing using phylogenies. *Bioinformatics.* 2005;21(5):676–679.
6. Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics.* 2009;25(11):1422–1423.
7. Tripp RA, et al. CX3C chemokine mimicry by respiratory syncytial virus G glycoprotein. *Nat Immunol.* 2001;2(8):732–738.
8. [RSV G central conserved domain structure / scFv 2D10; PDB 5WN9] — verify primary citation.
9. Wiley DC, Wilson IA, Skehel JJ. Structural identification of the antibody-binding sites of Hong Kong influenza haemagglutinin and their involvement in antigenic variation. *Nature.* 1981;289:373–378. (H3 antigenic sites)
10. The PyMOL Molecular Graphics System, open-source build, v3.1.0. Schrödinger, LLC.
11. Murrell B, et al. Detecting individual sites subject to episodic diversifying selection. *PLoS Genet.* 2012;8(7):e1002764. (MEME)

---

*Manuscript draft generated 2026-07-04 from the completed U6 analysis pipeline. Verify all references and fill placeholder author/affiliation/funding fields before submission.*
