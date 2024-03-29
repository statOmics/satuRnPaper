# satuRn: Scalable Analysis of Differential Transcript Usage for Bulk and Single-cell RNA-sequencing Applications

## Abstract

**Background**

Alternative splicing allows for the production of multiple, functionally distinct transcripts from a single genomic locus. However, the deregulation of this splicing process has been extensively reported as a cause of disease and as a hallmark of cancer. While several differential transcript usage (DTU) analysis tools already exist, they are not optimally suited to analyze large bulk RNA-seq and full-length scRNA-seq datasets, leaving the great potential of these data largely unexploited. 

**Results**

We therefore developed a novel DTU analysis tool named satuRn (acronym: Scalabale Analysis of differential Transcript Usage for RNa-seq data). 
We evaluated the performance of of satuRn on publicly available simulated and real bulk RNA-seq benchmark datasets, as well as on real scRNA-seq benchmark datasets. We found that the performance of satuRn was at least on par with the performances of the best tools from the literature. In addition, satuRn controls the FDR closer to the nominal level in general. 
We also show that satuRn scales towards the large data volumes generated by contemporary bulk and single-cell RNA-seq experiments, allowing for a transcriptome-wide analysis of datasets consisting of several thousands of cells in only a few minutes. Finally, we analyze a large full-length scRNA-seq case study dataset, where we obtain highly relevant biological results on isoform-level changes between cell types that would have remained obscured in a canonical differential gene expression (DGE) analysis.

**Conclusions**

We have shown that our novel DTU analysis tool, satuRn:
(i) is highly performant, 
(ii) provides a strict control of the FDR, 
(iii) scales seamlessly to the large data volumes of contemporary (sc-)RNA-seq datasets, 
(iv) allows for modelling complex experimental designs, 
(v) can deal with realistic proportions of zero counts and
(vi) provides direct inference on the biologically relevant transcript level.
To our knowledge, satuRn is the only DTU analysis method that combines all of the above qualities.

**Methods**

satuRn adopts a quasi-binomial generalized linear model framework. satuRn provides direct inference on DTU by modelling the relative usage of a transcript, in comparison to other transcripts from the same gene, between conditions of interest. To stabilize the estimation of the overdispersion parameter of the QB model, we borrow strength across transcripts by building upon the empirical Bayes methodology as introduced by Smyth et al. In order to control the number of false positive findings, an empirical null distribution is used to obtain the p-values, which are corrected for multiple testing with the FDR method of Benjamini and Hochberg. 

satuRn is implemented as an R software package and is available on Github from https://github.com/statOmics/satuRn.

***

## Availability of data

The datasets required to reproduce all results that are displayed in this publication (including supplementary materials) are available at Zenodo. This includes both the raw data and intermediate results. Note that at the top of each analysis script it is indicated which dataset is required as input for the script; it may thus not be necessary to download all datasets from [Zenodo](https://doi.org/10.5281/zenodo.4439415). For easily reproducing our analyses, place the downloaded data in the `Data` folder of your local clone of this repository.

***

## Analyses & Scripts

To reproduce the results that are displayed in this publication, proceed as follows:

1. Make a local clone of this Github repository
2. Open the selected R scripts (.Rmd) in the R project "DTU_paper.Rproj" in the root of this repository.
3. Look at the top of the R scripts what input data is required. Download this data from Zenodo (see Availability of data) and place the data in the `Data`  folder of this repository.
4. Run the analyses - the results will automatically be stored in the Git repository, either in the `Data` folder (if it is an intermediate results) or in the `Results` folder (if it is an end results, e.g. a figure). To improve the replicability of the analyses performed in this paper, a file `sessionInfo.txt` is included in the root of this repository, displaying the exact versions of each package and their respective dependencies that we have used for the analyses. 

This Github repository contains the following folders:

`CaseStudy`: This folder contains all the code for reproducing the case study analyses (analyses with satuRn, limma diffsplice and DoubleExpSeq).

`Data`: For reproducing the results displayed in this publication, input data for the different analysis should be downloaded from Zenodo (see Availability of data) and placed in the `Data` folder of **your local clone** of this Github page.

`Performance_benchmarks`: This folder contains all the code for reproducing the performance benchmark analyses.

`Results`:  This folder contains all the resulting figures, including those that were used in this publication and its supplementary materials.

`Scalability_benchmarks`: This folder contains all the code for reproducing the scalability benchmark analyses.

## Additional analysis in paper version 2 after reviewer comments

To accommodate the suggestions of the peer review, we included several 
additional analyses in the satuRn publication. Again, scripts are available from
Zenodo.

`CaseStudy_complex_design`: This folder contains a script for performing a data
analysis of the dataset of Tasic et al., but in contrast with the original case
study we here include the categorical covariate gender and the continuus 
covariate age in the specification of the DTU model.

`CaseStudy_EC`: Here, we again analyse the dataset by Tasic et al., but upon
quantifying the data on the level of equivalence classes as opposed to the 
transcript level.

`CaseStudy_exon`: Contains a script to analyse a dataset that has been 
previously quantified on the exon level.

`Scalability_benchmarks`: This folder has been extended to also perform a 
scalability benchmark of the different algorithms on bulk data.

`Sumstats`: This folder contains script to visualize and summarize some of the
key summary statistics of the different dataset, e.g., library size 
distributions, number of transcripts retained after filtering, etc. 


