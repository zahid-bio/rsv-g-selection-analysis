# Enrichment summary — U6 project

| Target | Region | n | hits | obs % | exp % | enrichment | p-value | sig |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| flu | H3 antigenic sites A-E | 5 | 5 | 100.0 | 39.43 | 2.54x | 0.0081 | YES |
| rsv_a (CCR) | RSV G CCR 157-198 | 3 | 0 | 0.0 | 13.9 | 0.0x | 1.0 | no |
| rsv_a (CX3C) | RSV G CX3C 182-186 | 3 | 0 | 0.0 | 1.58 | 0.0x | 1.0 | no |
| rsv_b (CCR) | RSV G CCR 157-198 | 3 | 0 | 0.0 | 13.9 | 0.0x | 1.0 | no |
| rsv_b (CX3C) | RSV G CX3C 182-186 | 3 | 0 | 0.0 | 1.58 | 0.0x | 1.0 | no |
| rsv_combined (CCR) | RSV G CCR 157-198 | 6 | 0 | 0.0 | 14.1 | 0.0x | 1.0 | no |
| rsv_combined (CX3C) | RSV G CX3C 182-186 | 6 | 0 | 0.0 | 1.66 | 0.0x | 1.0 | no |

All tests: 10,000 permutations, seed 42, one-sided (enrichment).

## Region-level selection (threshold-free companion)

Per-codon dN-dS (FEL beta-alpha) compared across regions using ALL codons (not just the FEL-significant handful). Negative effect size = CCR under stronger purifying selection.

| Target | Comparison | Test | Statistic | Effect size | p-value | sig |
|---|---|---|---|---|---:|:---:|
| rsv_a | Nterm/CCR/Mucin (omnibus) | Kruskal-Wallis | H=7.19 | eps2=0.018 | 0.0275 | YES |
| rsv_a | CCR < Mucin (more purifying) | Mann-Whitney U (1-sided) | U=1546 | rbc=-0.264 | 0.0066 | YES |
| rsv_a | CCR < non-CCR (mapped) | Mann-Whitney U (1-sided) | U=4172 | rbc=-0.224 | 0.0094 | YES |
| rsv_b | Nterm/CCR/Mucin (omnibus) | Kruskal-Wallis | H=2.76 | eps2=0.003 | 0.252 | no |
| rsv_b | CCR < Mucin (more purifying) | Mann-Whitney U (1-sided) | U=1705 | rbc=-0.154 | 0.0742 | no |
| rsv_b | CCR < non-CCR (mapped) | Mann-Whitney U (1-sided) | U=4494 | rbc=-0.151 | 0.0561 | no |

Rank-biserial (rbc): 0.1 small, 0.3 medium, 0.5 large. RSV-B A2 mapping loses 4 C-terminal mucin residues (a2_ref_len 294 vs 298; MAFFT --keeplength).
