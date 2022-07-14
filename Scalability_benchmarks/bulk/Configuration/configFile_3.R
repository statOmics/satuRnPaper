suppressPackageStartupMessages(library(rjson))

L <- list(countData = "./data/quantsf_counts_bulk.Rds", 
          tx2gene = "./data/tx2gene_bulk.Rds",
          resfilebase = "results",
          groupSizes = c(8,16,32), ## start low
          transcripts_amount = c(10000,30000))

countData_init <- readRDS(file = L$countData) # read in only once
tx2gene_init <- readRDS(file = L$tx2gene)

write(toJSON(L), file = "./Configuration/configFile.json")
print(L)
invisible(gc())


