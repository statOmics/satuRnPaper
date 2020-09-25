suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(limma))

run_edgeRDiffsplice <- function(L,countData,tx2gene) {
  message("edgeRDiffsplice")
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
    
    if(dim(quantsf_counts)[1] > i){
    
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
        	design <- sampleData[-which(colSums(quantsf_counts) == 0),]
        	quantsf_counts <- quantsf_counts[,-which(colSums(quantsf_counts) == 0)]
      	}

	y <- DGEList(counts = quantsf_counts, group = design$group)
        y <- calcNormFactors(y)
      
        design <- model.matrix(~group, data=design)
      
        y <- estimateDisp(y, design = design)
        fit <- glmFit(y, design = design) ## glmQLfit --> very poor results
      
        dtuEdgeR <- edgeR::diffSpliceDGE(
            glmfit = fit,
            exonid = tx2gene$TXNAME,
            geneid = tx2gene$GENEID,
            coef = "groupB",
            verbose = FALSE
        )
      
        localRes <- topSpliceDGE(
            dtuEdgeR,
            test = 'exon',
            number = Inf
        )

	finalResult <- nrow(localRes[which(localRes$FDR < 0.05),])
	
	mem <- gc()
    
    } else {
      ## make sure to have an output for result and info
        quantsf_counts <- matrix(data=NA, nrow=1,ncol=1)
        finalResult <- "More TXs asked than there were present"
	
	mem <- gc()
    }
  })
  
  output_limmaDiffsplice <- list(session_info = session_info,
                                 info = list(size = j, 
                                             TXs = nrow(quantsf_counts)),
                                 timing = timing,
				 removed = removed,
				 memory = mem,
                                 result = finalResult)
  
  assign("result", output_limmaDiffsplice, envir=globalenv())
  rm(list = ls())
  gc()
}
