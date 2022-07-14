# Scalability_benchmark folder

- To perform the scalability analyses, the unix executable `runtime` can 
be run from the command line using *./runtime*. There is a separte `runtime`
executable for the singlecell and bulk benchmarks, which can be found in the
`singlecell` and `bulk` subfolders, respectively. Running this script requires 
having the `Scalability_data` folder, which can be downloaded from Zenodo, 
to be in the `Data` folder of the repository (unzipped). `runtime` performs the
analyses on all scalability benchmark data for all DTU methods assessed in this 
publication. In total, the scalability benchmark ran for approximately 29 hours 
on our system for the single-cell analysis.
- Intermediate results (RData files containing the runtimes) are available from 
the `Results < Scalability_benchmarks` folder, subfolders `singlecell` and 
`bulk`.
- To generate the figures with respect to scalability displayed in this
benchmark, run the `scalabilityBenchmark.Rmd` script.

Note: running BANDITS requires equivalence class count files, which are 27Gb in 
size. We did not upload these files to zenodo, but they can be requested.