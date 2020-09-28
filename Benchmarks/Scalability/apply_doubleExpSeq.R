suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DoubleExpSeq))

run_doubleExpSeq <- function(L,countData,tx2gene) {
  message("doubleExpSeq")
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

	tx2gene <- tx2gene[which(tx2gene$TXNAME %in% TXnames),]
    
        design <- as.data.frame(matrix(nrow = ncol(quantsf_counts), ncol=2))
        colnames(design) <- c("sampleId", "group")
        design$sampleId <- colnames(quantsf_counts)
        design$group <- sample(rep(c("A","B"),each=i))
	
	## check if there are cells without counts due to filtering/subsampling
	## if so (unlikely); remove them, but store the fact they are removed
	if(length(which(colSums(quantsf_counts) == 0)) > 0){
        	removed <- length(which(colSums(quantsf_counts) == 0))
        	design <- design[-which(colSums(quantsf_counts) == 0),]
        	quantsf_counts <- quantsf_counts[,-which(colSums(quantsf_counts) == 0)]
      	}

	# get tx2gene in better format
        geneForEachTx <- tx2gene$GENEID[match(rownames(quantsf_counts),tx2gene$TXNAME)]
        groupID <- as.character(geneForEachTx)
      
        ## get the "total" counts
        stopifnot(class(quantsf_counts) %in% c("matrix", "data.frame"))
        countData <- as.matrix(quantsf_counts)
        stopifnot(class(groupID) %in% c("character", "factor"))
        stopifnot(length(groupID) == nrow(quantsf_counts))

        forCycle <- split(1:nrow(quantsf_counts), as.character(groupID))
        all <- lapply(forCycle, function(i) {
            sct <- quantsf_counts[i, , drop = FALSE]
            rs <- t(sapply(1:nrow(sct), function(r) colSums(sct[, , drop = FALSE])))
            rownames(rs) <- rownames(sct)
            rs
        })
        totalCount <- do.call(rbind, all)
        stopifnot(all(rownames(quantsf_counts) %in% rownames(totalCount)))
        totalCount <- totalCount[rownames(quantsf_counts), ]
        
        group <- as.factor(design$group)
        
        mDouble <- suppressWarnings(DBGLM1(y = quantsf_counts,
                  m = totalCount,
                  groups = group))

	finalResult <- nrow(mDouble$Sig[which(mDouble$Sig[,4] < 0.05),])
	
	mem <- gc()
    
    } else {
      ## make sure to have an output for result and info
        quantsf_counts <- matrix(data=NA, nrow=1,ncol=1)
        finalResult <- "More TXs asked than there were present"
	
	mem <- gc()
    }
  })
  
  output_DoubleExpSeq <- list(session_info = session_info,
                                 info = list(size = i, 
                                             TXs = nrow(quantsf_counts)),
                                 timing = timing,
				 removed = removed,
				 memory = mem,
                                 result = finalResult)
  
  assign("result", output_DoubleExpSeq, envir=globalenv())
  rm(list = ls())
  gc()
}
