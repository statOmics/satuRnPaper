suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DRIMSeq))

run_DRIMSeq <- function(L,countData,tx2gene) {
  message("DRIMSeq")
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
        dim(quantsf_counts)
    
        design <- as.data.frame(matrix(nrow = ncol(quantsf_counts), ncol=2))
        colnames(design) <- c("sample_id", "group")
        design$sample_id <- colnames(quantsf_counts)
        design$group <- sample(rep(c("A","B"),each=L$groupSizes[i])) ## shouldn't this be i in stead of L$groupSizes[i]?
    
        geneForEachTx <- tx2gene[match(rownames(quantsf_counts),tx2gene$TXNAME),2]
	
	## check if there are cells without counts due to filtering/subsampling
	## if so (unlikely); remove them, but store the fact they are removed
	if(length(which(colSums(quantsf_counts) == 0)) > 0){
        	removed <- length(which(colSums(quantsf_counts) == 0))
        	design <- sampleData[-which(colSums(quantsf_counts) == 0),]
        	quantsf_counts <- quantsf_counts[,-which(colSums(quantsf_counts) == 0)]
      	}
    
        ## To make DRIMSeq analysis reproducible
        set.seed(124)
      
        # get tx2gene in better format
        quantsf_counts <- round(quantsf_counts)
        quantsf_counts <- as.data.frame(quantsf_counts)

        geneForEachTx <- tx2gene$GENEID[match(rownames(quantsf_counts),tx2gene$TXNAME)]

        quantsf_counts$gene_id <- geneForEachTx
        quantsf_counts$feature_id <- row.names(quantsf_counts)
        rownames(quantsf_counts) <- c()

        d <- dmDSdata(counts = quantsf_counts, samples = design)

        design <- model.matrix(~ group, data = DRIMSeq::samples(d))

        d <- dmPrecision(d, design = design, verbose=0)
        d <- dmFit(d, design = design, verbose=0)
        d <- dmTest(d, coef = 2, verbose=0)

        localRes <- DRIMSeq::results(d, level = "feature")[, c(
            'feature_id',
            'gene_id',
            'lr',
            'df',
            'pvalue',
            'adj_pvalue'
        )]
	localRes <- as.data.frame(localRes)
      
        finalResult <- nrow(localRes[which(localRes$adj_pvalue < 0.05),])
	
	mem <- gc()
    
    } else {
      ## make sure to have an output for result and info
        quantsf_counts <- matrix(data=NA, nrow=1,ncol=1)
        finalResult <- "More TXs asked than there were present"
	
	mem <- gc()
    }
  })
  
  output_DRIMSeq <- list(session_info = session_info,
                                 info = list(size = i, 
                                             TXs = nrow(quantsf_counts)),
                                 timing = timing,
				 removed = removed,
				 memory = mem,
                                 result = finalResult)
  
  assign("result", output_DRIMSeq, envir=globalenv())
  rm(list = ls())
  gc()
}
