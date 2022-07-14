suppressPackageStartupMessages(library(rjson))

L <- list(countData = "./Data/Scalability_data/quantsf_counts.Rds", 
          tx2gene = "./Data/Scalability_data/tx2gene.Rds",
          resfilebase = "results",
          groupSizes = c(16),
          transcripts_amount = c(2000, 5000, 15000, 20000, 25000, 35000))

countData_init <- readRDS(file = L$countData) # read in only once
tx2gene_init <- readRDS(file = L$tx2gene)

write(toJSON(L), file = "./Scalability_benchmark/Configuration/configFile.json")
print(L)
invisible(gc())


