#!/usr/bin/env bash
# Recombination screen (GARD, Genetic Algorithm for Recombination Detection)
# for the RSV-A and RSV-B G codon alignments.
#
# GARD requires a working MPI environment. Two supported routes:
#
# (1) Datamonkey web server (recommended, no local setup):
#       https://www.datamonkey.org/gard
#     Upload results/rsv_a_codon_aligned.fasta and results/rsv_b_codon_aligned.fasta,
#     genetic code = Universal, site-to-site rate variation = None/GDD as offered,
#     and record: number of inferred breakpoints, their locations, and whether the
#     resulting partition topologies are incongruent (KH test) at p < 0.05.
#
# (2) Local HyPhy-MPI (requires a functioning MPI stack; the biocontainers image
#     used elsewhere in this project did not provide a working MPI communicator):
#
#   for k in rsv_a rsv_b; do
#     mpirun -np 4 HYPHYMPI gard \
#       --alignment results/${k}_codon_aligned.fasta \
#       --output    results/${k}_gard.json
#   done
#
# Report in Methods/Results: number of breakpoints and whether any breakpoint is
# associated with a statistically supported change in tree topology (which would
# motivate partitioned selection analyses). If no topology-incongruent breakpoint
# is found, the single-tree FEL/MEME analyses are appropriate as run.
set -euo pipefail
echo "See comments in this script. Preferred route: Datamonkey (https://www.datamonkey.org/gard)."
