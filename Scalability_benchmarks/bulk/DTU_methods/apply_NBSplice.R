suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(NBSplice))

run_NBSplice <- function(L,countData,tx2gene) {
  message("NBSplice")
  session_info <- sessionInfo()
  timing <- system.time({
    
    tx2gene <- tx2gene_init
    quantsf_counts <- countData_init
    removed <- "none"
    
    ## select cells
    set.seed(123)
    selected_cells <- sample(colnames(quantsf_counts), 2*i)
    head(selected_cells)

    quantsf_counts <- quantsf_counts[,selected_cells]
    
    ## filter TXs with zero rowSum after cell selection
    quantsf_counts <- quantsf_counts[which(rowSums(quantsf_counts)!=0),]
    dim(quantsf_counts)
    tx2gene <- tx2gene[tx2gene$TXNAME %in% rownames(quantsf_counts),]
    
    if(dim(quantsf_counts)[1] > j){
    
        ## select transcripts
        k <- 1
        randomized <- sample(unique(tx2gene$GENEID))
        TXnames <- c()

        while (length(TXnames) < j){
            TXnames <- c(TXnames, tx2gene$TXNAME[tx2gene$GENEID == randomized[k]])
            k <- k + 1
        }

        quantsf_counts <- quantsf_counts[TXnames,]
        quantsf_counts <- as.data.frame(quantsf_counts)
    
        sampleData <- as.data.frame(matrix(nrow = ncol(quantsf_counts), ncol=2))
        colnames(sampleData) <- c("id", "condition")
        sampleData$id <- colnames(quantsf_counts)
        sampleData$condition <- sample(rep(c("A","B"),each=i))
	sampleData$condition <- as.factor(sampleData$condition)
	rownames(sampleData) <- sampleData$id
    
        geneForEachTx <- tx2gene[match(rownames(quantsf_counts),tx2gene$TXNAME),2]
	
	## check if there are cells without counts due to filtering/subsampling
	## if so (unlikely); remove them, but store the fact they are removed
	if(length(which(colSums(quantsf_counts) == 0)) > 0){
        	removed <- length(which(colSums(quantsf_counts) == 0))
        	design <- sampleData[-which(colSums(quantsf_counts) == 0),]
        	quantsf_counts <- quantsf_counts[,-which(colSums(quantsf_counts) == 0)]
      	}
      
        # get tx2gene in better format
        quantsf_counts <- round(quantsf_counts)
        quantsf_counts <- as.data.frame(quantsf_counts)

	geneIso <- as.data.frame(tx2gene[,c(2,1)])
	rownames(geneIso) <- geneIso[,2]
	colnames(geneIso) <- c("gene_id", "isoform_id")
	geneIso <- geneIso[which(geneIso$isoform_id %in% TXnames),]
	geneIso <- geneIso[match(rownames(quantsf_counts), geneIso$isoform_id),]

	colName <- "condition"

	quantsf_counts <- round(quantsf_counts)

	myIsoDataSet <- IsoDataSet(quantsf_counts, sampleData, colName, geneIso)
	myDSResults <- NBTest(object=myIsoDataSet, colName=colName, test="F")

	localRes <- myDSResults@results
	localRes <- localRes[,c(1,2,7)]
	colnames(localRes) <- c("TXNAME", "GENEID", "p_value")
	localRes <- as.data.frame(localRes)
	localRes$adj_pvalue <- p.adjust(localRes$p_value, method="BH")

        finalResult <- nrow(localRes[which(localRes$adj_pvalue < 0.05),])
	
	mem <- gc()
    
    } else {
      ## make sure to have an output for result and info
        quantsf_counts <- matrix(data=NA, nrow=1,ncol=1)
        finalResult <- "More TXs asked than there were present"
	
	mem <- gc()
    }
  })
  
  output_NBSplice <- list(session_info = session_info,
                                 info = list(size = i, 
                                             TXs = nrow(quantsf_counts)),
                                 timing = timing,
				 removed = removed,
				 memory = mem,
                                 result = finalResult)
  
  assign("result", output_NBSplice, envir=globalenv())
  rm(list = ls())
  gc()
}
