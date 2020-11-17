suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(locfdr))
suppressPackageStartupMessages(library(satuRn))

## run satuRn
run_satuRn <- function(L,countData,tx2gene) {
  message("satuRn")
  session_info <- sessionInfo()
  timing <- system.time({

    tx2gene <- tx2gene_init
    quantsf_counts <- countData_init
    removed <- "none"
    
    ## select cells
    set.seed(123)
    selected_cells <- sample(colnames(quantsf_counts), 2*i)

    quantsf_counts <- quantsf_counts[,selected_cells]
    
    ## filter TXs with zero rowSum after cell selection
    quantsf_counts <- quantsf_counts[which(rowSums(quantsf_counts)!=0),]
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
    
        design <- as.data.frame(matrix(nrow = ncol(quantsf_counts), ncol=2))
        colnames(design) <- c("sampleId", "group")
        design$sampleId <- colnames(quantsf_counts)
        #design$group <- sample(rep(c("A","B"),each=L$groupSizes[i]))
        design$group <- sample(rep(c("A","B"),each=i))
	
	      ## check if there are cells without counts due to filtering/subsampling
	      ## if so (unlikely); remove them, but store the fact they are removed as it may impact speed
	      if(length(which(colSums(quantsf_counts) == 0)) > 0){
        	  removed <- length(which(colSums(quantsf_counts) == 0))
        	  design <- design[-which(colSums(quantsf_counts) == 0),]
        	  quantsf_counts <- quantsf_counts[,-which(colSums(quantsf_counts) == 0)]
      	}

	      colnames(tx2gene)[1:2] <- c("isoform_id", "gene_id")
        group <- factor(design$group)
        tx2gene <- tx2gene[which(tx2gene$isoform_id %in% rownames(quantsf_counts)),]
        
        sumExp <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=quantsf_counts), colData = design, rowData = tx2gene)
        MultiAssayExperiment::metadata(sumExp)$formula <- ~ 0 + group
      
        
        sumExp <- satuRn:::fitQB(object = sumExp,
                                 parallel = FALSE,
                                 BPPARAM = BiocParallel::bpparam(),
                                 verbose = TRUE
        )
        
        L <- matrix(0, ncol = 1, nrow = ncol(design))
        rownames(L) <- colnames(design)
        colnames(L) <- c("AvsB")
        L[1:2,] <- c(-1,1)
        
        sumExp <- satuRn::testDTU(object = sumExp, contrasts = L[,1,drop=FALSE], plot=F, sort = F)
        
        output <- data.frame(rownames(quantsf_counts))
        rownames(output) <- rownames(quantsf_counts)
        colnames(output) <- "TXNAME"
        
        output$GENEID <- tx2gene[match(output$TXNAME,tx2gene$TXNAME),"GENEID"]
        output$tvalue <- rowData(sumExp)[["fitQBResult_AvsB"]]$t
        output$p_value_raw <- rowData(sumExp)[["fitQBResult_AvsB"]]$pval
        output$FDR_raw <- rowData(sumExp)[["fitQBResult_AvsB"]]$regular_FDR
        output$p_value <- rowData(sumExp)[["fitQBResult_AvsB"]]$empirical_pval
        output$FDR <- rowData(sumExp)[["fitQBResult_AvsB"]]$empirical_FDR
    
        finalResult <- nrow(output[which(output$FDR < 0.05),])
	
	      mem <- gc()
    } else {
      ## make sure to have an output for result and info
        quantsf_counts <- matrix(data=NA, nrow=1,ncol=1)
        finalResult <- "More TXs asked than there were present"
	      mem <- gc()
    }
  })
  
  output_satuRn <- list(session_info = session_info,
                                 info = list(size = i, 
                                             TXs = nrow(quantsf_counts)),
                                 timing = timing,
				                         removed = removed,
				                         memory = mem,
                                 result = finalResult)
  
  assign("result", output_satuRn, envir=globalenv())
  rm(list = ls())
  gc()
}
