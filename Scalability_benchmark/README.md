## Scalability_benchmark folder

## Current problem: BANDITS requires ECC files, which are the bulk of the 27Gb of the salmon files --> no solution, so either ECCs must be provided in Zenodo or the BANDITS analysis cannot be reproduced. I tink we should upload to Zenodo.

- To perform the scalability analyses, the unix executable `runtime` can simply be run from the command line using ./runtime. Running this script requires having the `Scalability_runtimes` folder, which can be downloaded from Zenodo, to be in the `Data` folder of the repositroy (unzipped).  `runtime` performs the analyses on all scalability benchmark data for all DTU methods assessed in this publication. In total, the scalability benchmark ran for approximately 29 hours on our system.
- To generate the figures with respect to scalability displayed in this benchmark, run the `scalabilityBenchmark.Rmd` script.

