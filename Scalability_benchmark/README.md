# DTU_paper: scalability benchmark

The script initialFiltering.R is the first to run. It generates the required count matrix (quantsf_counts.Rds),
a file linking transcripts to genes (tx2gene.Rds) and the estimated effective transcripts lenghts (eff_len.Rds).
The latter file is required to perform a BANDITS analysis.

The file