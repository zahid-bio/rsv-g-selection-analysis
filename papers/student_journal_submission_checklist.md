# Student Journal Submission Checklist

Manuscript: **RSV G selection analyses support constraint of the CX3C/cystine-noose core and diversification in mucin-like regions**
Author: Zahid Ayomide Nassoro-Ally (Independent Researcher, Chicago, IL, USA)
ORCID: [0009-0002-0550-5115](https://orcid.org/0009-0002-0550-5115)
Correspondence: nassorozahid@gmail.com

Status legend: [x] done · [ ] to do by author · [~] partial / needs external step

## Manuscript package

- [x] Manuscript source present (`papers/RSV_G_SELECTION_PAPER.tex`).
- [~] Final PDF built successfully (`papers/RSV_G_Selection_Paper.pdf`) — rebuild after any late edits; see "Building the PDF" below.
- [x] Title in PDF matches repository README title.
- [x] Author block correct: name, affiliation, ORCID, correspondence email.
- [x] Abstract present and self-contained.
- [x] Keywords present.
- [x] Running title present.

## Figures and tables

- [x] All figures referenced in text and present in `figures/` (Figures 1–4; Supplementary S1–S5).
- [x] All main tables present (Table 1 enrichment; Table 2 FEL/MEME region counts; Table 3 region-level dN−dS).
- [x] Supplementary tables labeled S1, S2 (counter reset, "Supplementary Table" prefix).
- [x] Every figure/table cited in order in the text.
- [x] Captions understandable without reading the full paper.
- [x] Figure 1 (flu HA validation) shows p = 0.0081 with observed/expected labels.
- [x] Figure 2 (RSV lollipop + CCR/CX3C shading) present.
- [x] Figure 3 (region-level dN−dS, RSV-A and RSV-B) present.
- [x] Figure 4 (PDB 5WN9 structure panel) present.
- [ ] If journal requires separate high-resolution figure files, export from `figures/` (PNG present; regenerate as PDF/TIFF if required).

## Numbers and reproducibility

- [x] Sequence counts match files: RSV-A n = 108, RSV-B n = 115, flu H3 HA n = 148.
- [x] Tool versions in Methods match README (MAFFT 7.526; IQ-TREE 3.1.3; HyPhy 2.5.100 Docker; Python 3.11; Biopython 1.87; PyMOL 3.1.0).
- [x] ModelFinder models match (TVM+F+I+R2 flu; TN+F+I+R2 RSV-A and RSV-B).
- [x] Region coordinates correct (CCR 157–198; CX3C 182–186; N-term 1–156; mucin 199–298; 5WN9 ordered core 169–189).
- [x] Permutation test described (one-sided, 10,000 permutations, seed 42).
- [x] FEL/MEME threshold stated as unadjusted p ≤ 0.05 and interpreted cautiously.
- [x] Every table value matches its source CSV in `results/` (verified programmatically).
- [x] SHA256 checksums recorded (`notes/lab_notebook.txt`; `results/meme/SHA256SUMS`).

## References

- [x] All in-text citations have a bibliography entry (16/16).
- [x] All bibliography entries are cited (16/16).
- [x] PDB 5WN9 cited (Fedechkin et al. 2018 + PDB DOI).
- [x] PyMOL cited.
- [x] Influenza HA antigenic-site definitions cited (Wiley/Wilson/Skehel 1981).
- [x] DOIs formatted consistently.
- [ ] Reference list spot-checked against primary sources by author before submission.

## Honesty and scope

- [x] No overclaiming language ("prove", "demonstrate definitively", "validated target" absent).
- [x] Conservative verbs used ("support", "consistent with", "suggest", "do not by themselves establish").
- [x] Limitations section complete (power, subtype asymmetry, mapping loss, structural coverage, post-hoc localization, sampling metadata, multiplicity, distributional assumptions, recombination/tree uncertainty, method scope).
- [x] No TODO/DRAFT wording in Data and code availability.
- [x] Conflicts of interest stated (none).
- [x] Funding stated (none).
- [x] Acknowledgements stated (none — no external help to acknowledge).

## Licensing and repository

- [x] Code license present (MIT, `LICENSE_CODE`).
- [x] Manuscript/data license present (CC BY 4.0, `LICENSE_DATA_MANUSCRIPT`).
- [x] Third-party PDB 5WN9 coordinates not claimed as own work.
- [x] README license section finalized (no "to be added").
- [ ] Repository made public on GitHub before/at submission.
- [ ] Tagged GitHub release created (e.g., v0.1-preprint) so the submitted state is frozen.
- [ ] Optional: archive the release on Zenodo and add the DOI to the manuscript.

## Author actions before clicking submit

1. Rebuild the PDF from the final `.tex` (see below) and confirm no unresolved references/citations.
2. Confirm correspondence email `nassorozahid@gmail.com` is monitored.
3. Decide whether the target student journal requires an AI-use statement; if so, paste the text from `papers/AI_use_statement.md` into the manuscript (see task note there).
4. Make the GitHub repository public and create a tagged release.
5. Upload manuscript PDF + any separately required figure/supplementary files to the journal portal.
6. Fill in the journal name and editor name in `papers/student_journal_cover_letter.md`.

## Building the PDF

The manuscript is plain LaTeX (`article` class) with a `thebibliography` block (no external `.bib`, no bibtex run required). Two build options:

Local (if a TeX distribution is installed):

```
latexmk -pdf -outdir=papers papers/RSV_G_SELECTION_PAPER.tex
```

Docker (no local TeX install needed):

```
docker run --rm -v "$PWD":/work -w /work texlive/texlive:latest \
  latexmk -pdf -interaction=nonstopmode -outdir=papers papers/RSV_G_SELECTION_PAPER.tex
```

After building, confirm the log shows no `LaTeX Warning: Reference ... undefined` and no `Citation ... undefined`.
