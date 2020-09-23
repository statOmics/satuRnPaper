dir <- "/home/compomics/salmon_Chen_new/"
samples <- read.table(file.path(dir, "SRRIDChen.txt"), header = FALSE)

files <- file.path(dir, samples$V1, "quant.sf")
names(files) <- paste0("sample", 1:length(samples$V1))
all(file.exists(files))

library(AnnotationHub)
## Load the annotation resource.
ah <- AnnotationHub()

all <- query(ah, "EnsDb")

library(ensembldb)
ahEdb <- all[["AH75036"]] #for Mus musculus

## retrieve all txs (and corresponding genes)
txs <- transcripts(ahEdb)

rm(ah,all)
gc()

tx2gene <- as.data.frame(matrix(data = NA, nrow = length(txs), ncol = 2))
colnames(tx2gene) <- c("TXNAME","GENEID")
tx2gene$TXNAME <- txs$tx_id
tx2gene$GENEID <- txs$gene_id

library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = TRUE)

quantsf_counts <- txi$counts
sum(quantsf_counts==0)/(sum(quantsf_counts==0)+sum(quantsf_counts!=0))
dim(quantsf_counts)

eff_len <- BANDITS::eff_len_compute(x_eff_len = txi$length) # save for BANDITS

gc()

metaData <- read.csv("/home/compomics/scTime_Chen/Chen_metadata.csv")
selected <- unname(substring(files,34,43))
colnames(quantsf_counts) <- selected

table(interaction(metaData$Cell_type, metaData$sex))

tx2gene <- tx2gene[which(tx2gene$TXNAME %in% sub("\\..*", "", rownames(quantsf_counts))),]

tx2gene$TXNAME <- rownames(quantsf_counts)[match(tx2gene$TXNAME,sub("\\..*", "", rownames(quantsf_counts)))]

#Filter transcripts and cells, so that the transcripts that will be selected later have at least a minimum in expression and the cells have a certain minimal quality. Otherwise, the scalability profiling may not be representative, i.e. if  the poor quality of cells/transcripts (that in real analyses would be filtered out) are causing artificial problems in parameter estimation.

#1. Filter out genes with a very low count (very lenient filter)

library(edgeR)
filter <- filterByExpr(quantsf_counts, min.count = 1, min.total.count = 1, large.n = 1, min.prop = 0.1)
summary(filter)

quantsf_counts <- quantsf_counts[filter,]
dim(quantsf_counts)

tx2gene <- tx2gene[tx2gene$TXNAME %in% rownames(quantsf_counts),]

#2. Filter out genes with only one Tx

tx2gene <- tx2gene[which(tx2gene$GENEID %in% tx2gene$GENEID[duplicated(tx2gene$GENEID)]),]

quantsf_counts <- quantsf_counts[rownames(quantsf_counts) %in% tx2gene$TXNAME,]
dim(quantsf_counts)

#3. Filter out cells with low library size

#hist(colSums(quantsf_counts), breaks =100)
#abline(v=750000, col="red", lwd=2)

quantsf_counts <- quantsf_counts[,colSums(quantsf_counts) > 750000]
dim(quantsf_counts)

#4. Filter out cells with low transcriptome complexity (unique TXs)

unique_TXs_sample <- vector(mode = "numeric",length = dim(quantsf_counts)[2])

for (i in 1:dim(quantsf_counts)[2]){
  unique_TXs_sample[i] <- length(which(quantsf_counts[,i] > 0)) 
}

#hist(
 #   unique_TXs_sample,
  #  breaks = 100, xlab = "unique_TXs_per_cell"
#)
#abline(v = 12500, col = "red", lwd=2)

quantsf_counts <- quantsf_counts[,-which(unique_TXs_sample < 12500)]
dim(quantsf_counts)

#These filtering steps (very lenient) lead to a dataset of 551 cells and 37.389 transcripts.

saveRDS(quantsf_counts,"/home/compomics/scTime_Chen/quantsf_counts.Rds")
saveRDS(tx2gene,"/home/compomics/scTime_Chen/tx2gene.Rds")
saveRDS(eff_len, file = "/home/compomics/scTime_Chen/eff_len.Rds") # save for bandits




