suppressPackageStartupMessages(library(rjson))

L <- list(countData = "./data/quantsf_counts_bulk.Rds", 
          tx2gene = "./data/tx2gene_bulk.Rds",
          resfilebase = "results",
          groupSizes = c(16),
          transcripts_amount = c(2000, 5000, 15000, 20000, 25000, 35000))

countData_init <- readRDS(file = L$countData) # read in only once
tx2gene_init <- readRDS(file = L$tx2gene)

write(toJSON(L), file = "./Configuration/configFile.json")
print(L)
invisible(gc())


