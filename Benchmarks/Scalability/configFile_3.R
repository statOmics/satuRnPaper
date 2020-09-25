suppressPackageStartupMessages(library(rjson))

L <- list(countData = "/home/compomics/scTime_Chen/quantsf_counts.Rds", 
          tx2gene = "/home/compomics/scTime_Chen/tx2gene.Rds",
          resfilebase = "results",
          groupSizes = c(8,16,32), ## start low
          transcripts_amount = c(10000,30000))

countData_init <- readRDS(file = L$countData) # read in only once
tx2gene_init <- readRDS(file = L$tx2gene)

write(toJSON(L), file = "config/configFile.json")
print(L)
gc()


