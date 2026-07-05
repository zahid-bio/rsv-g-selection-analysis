# GitHub Release Checklist for v0.1-preprint

**Status:** Ready for release  
**Current commits ahead of origin/main:** 2  
**Latest commit:** 039d912 (Add publication-ready licenses and finalize manuscript metadata)

---

## Pre-Release Verification ✅

- ✅ **Licenses added:**
  - `LICENSE_CODE` — MIT License for code (code/*.py)
  - `LICENSE_DATA_MANUSCRIPT` — CC-BY 4.0 for manuscript/figures/data

- ✅ **Manuscript finalized:**
  - Data and code availability: license terms clarified
  - Reproducibility note: DRAFT status removed
  - Acknowledgements: "None. All analyses performed independently."
  - Email: zahidtzzzz@gmail.com confirmed

- ✅ **Documentation:**
  - README title: matches manuscript title
  - Reproducibility notes: MEME_REPRODUCIBILITY_QUICKSTART.md, full checklist, summary
  - SHA256SUMS: all outputs checksummed (results/meme/SHA256SUMS)

- ✅ **Git commits:**
  - 5a5f78a: Make MEME analysis fully reproducible
  - 039d912: Add publication-ready licenses

---

## Step 1: Push Commits to GitHub

```bash
git push origin main
```

Expected output:
```
To https://github.com/zahid-bio/rsv-g-selection-analysis.git
   [hash] main -> main
```

---

## Step 2: Create GitHub Release

### Via GitHub Web UI:

1. Go to: https://github.com/zahid-bio/rsv-g-selection-analysis
2. Click **"Releases"** (right sidebar)
3. Click **"Create a new release"**
4. Fill in the form:

   **Tag version:** `v0.1-preprint`
   
   **Release title:** 
   ```
   RSV G selection analyses v0.1-preprint
   ```
   
   **Description:**
   ```
   Pre-publication release of the RSV G selection analysis manuscript and pipeline.
   
   ## What's Included
   
   - Manuscript: "RSV G selection analyses support CX3C/cystine-noose constraint and diversification in mucin-like regions"
   - Reproducible analysis pipeline: FEL (Fixed Effects Likelihood) and MEME (Mixed Effects Model of Evolution)
   - Raw outputs: MEME JSON files with SHA256 checksums
   - Parsed summaries: Regional site counts, FEL/MEME comparison tables
   - Complete documentation: Reproducibility notes, verification checklist, quick-start guide
   
   ## Licenses
   
   - Code (Python scripts): MIT License
   - Manuscript, figures, and data: CC-BY 4.0
   
   ## Inputs
   
   - Influenza H3N2 HA: 148 sequences, 578 codons
   - RSV-A G: 108 sequences, 322 codons
   - RSV-B G: 115 sequences, 336 codons
   
   All alignments codon-aware (MAFFT v7.526); phylogenies inferred with IQ-TREE v3.1.3.
   FEL and MEME analyses with HyPhy v2.5.100 (Docker image: quay.io/biocontainers/hyphy:2.5.100--h74d3ee0_0).
   
   ## Key Findings
   
   - CX3C motif is constrained: 0 FEL-positive sites, 0 MEME-positive sites in either RSV subtype
   - Broader CCR has sparse episodic signals: 1 MEME-positive site per subtype (non-CX3C)
   - Mucin-like domain is selection hotspot: 4-5 MEME sites per subtype
   - Influenza control: 6/11 MEME sites in antigenic regions (validates pipeline)
   
   ## Reproducibility
   
   All outputs are checksummed and fully reproducible from frozen inputs:
   - See MEME_REPRODUCIBILITY_QUICKSTART.md for 3-step verification
   - See papers/reproducibility_note.md for complete documentation
   - See results/meme/SHA256SUMS for checksum verification
   
   ## Citation
   
   Zahid. "RSV G selection analyses support CX3C/cystine-noose constraint and diversification in mucin-like regions." 
   GitHub repository: https://github.com/zahid-bio/rsv-g-selection-analysis (v0.1-preprint).
   
   ## Next Steps
   
   - Submit to bioRxiv or medRxiv as preprint
   - Archive on Zenodo after preprint posting (will add DOI to manuscript)
   - Update this release with Zenodo link after archival
   ```

5. Choose: **"This is a pre-release"** (checkbox) ✅
6. Click **"Publish release"**

### Via GitHub CLI (if installed):

```bash
gh release create v0.1-preprint \
  --title "RSV G selection analyses v0.1-preprint" \
  --notes-file RELEASE_NOTES.md \
  --prerelease
```

---

## Step 3: Verify Release

After creating the release, verify:

1. ✅ Visit: https://github.com/zahid-bio/rsv-g-selection-analysis/releases
2. ✅ Confirm v0.1-preprint appears with correct title and description
3. ✅ Confirm tag points to commit 039d912

---

## Step 4: Add Zenodo After Preprint Posting

After posting to bioRxiv/medRxiv:

1. Go to https://zenodo.org and log in (or create account if needed)
2. Create new deposit with:
   - Title: "RSV G selection analyses support CX3C/cystine-noose constraint and diversification in mucin-like regions"
   - Description: Link to preprint and GitHub release
   - Upload: Zip file of this repository (or just key files)
3. Get Zenodo DOI and update manuscript:

```markdown
## Data and code availability

All sequence accessions, frozen FASTA inputs, alignments, trees, FEL outputs, MEME outputs, analysis scripts, 
figures, and the lab notebook with SHA256 checksums are provided in the project repository: 
https://github.com/zahid-bio/rsv-g-selection-analysis (Zenodo DOI: 10.5281/zenodo.XXXXXXX).
```

4. Also update this release description with Zenodo link

---

## Commit Summary

```
5a5f78a - Make MEME analysis fully reproducible with comprehensive documentation
  - SHA256SUMS for all MEME outputs
  - MEME reproducibility checklist and summary
  - Updated papers/reproducibility_note.md with HyPhy versions

039d912 - Add publication-ready licenses and finalize manuscript metadata
  - LICENSE_CODE (MIT) for Python scripts
  - LICENSE_DATA_MANUSCRIPT (CC-BY 4.0) for manuscript/data
  - Updated manuscript Data and code availability section
  - Updated README title to match manuscript title
```

---

## Release Notes Template

**For bioRxiv/medRxiv submission:**

```
Title: RSV G selection analyses support CX3C/cystine-noose constraint and diversification in mucin-like regions

Authors: Zahid

Data and Code Availability:
All analysis code, data, and reproducibility documentation are available at:
- GitHub: https://github.com/zahid-bio/rsv-g-selection-analysis (v0.1-preprint)
- Zenodo: [DOI to be added after release]

Licenses:
- Code: MIT License
- Manuscript, figures, and data: CC-BY 4.0

Key Reproducibility Features:
- All MEME outputs checksummed (SHA256) in results/meme/SHA256SUMS
- Parsing script regenerates all CSVs byte-for-byte from raw JSON
- Complete HyPhy version specification and command lines included
- See MEME_REPRODUCIBILITY_QUICKSTART.md for 3-step verification
```

---

## Post-Release Checklist

- [ ] Verify release is visible at https://github.com/zahid-bio/rsv-g-selection-analysis/releases
- [ ] Update manuscript with Zenodo DOI (after archival)
- [ ] Add link to this release in acknowledgements or data availability section
- [ ] Share release link with collaborators
- [ ] Submit to preprint server with GitHub/Zenodo DOI
- [ ] Monitor release for issues; create patch releases if needed

---

## Release v0.1-preprint Status

**Ready to publish immediately.** All preparation complete:

✅ Code licensed (MIT)  
✅ Data/manuscript licensed (CC-BY 4.0)  
✅ Reproducibility documented and verified  
✅ Manuscript finalized (title, email, licenses)  
✅ Git commits prepared  
✅ Ready for GitHub release and preprint posting  

**Next action:** `git push origin main` and create release on GitHub.
