# DTU_paper: scalability benchmark

Current problems;
1. initialFiltering.R requires the salmon files (27Gb) --> can be resolved when only providing count matrix after tximport
2. BANDITS requires ECC files, which are the bulk of the 27Gb of the salmon files --> no solution

The script initialFiltering.R is the first to run. It generates the required count matrix (quantsf_counts.Rds),
a file linking transcripts to genes (tx2gene.Rds) and the estimated effective transcripts lenghts (eff_len.Rds).
The latter file is required to perform a BANDITS analysis.

The file