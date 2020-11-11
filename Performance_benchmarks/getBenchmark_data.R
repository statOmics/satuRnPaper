getBenchmark_data <- function(countData, metaData, filter,edgeR_filter_spec, nrRepList, fracGenesAffected){
  
  lapply(c(1:length(nrRepList)), function(x) {
    
    # Step 1: Extract random sub-sample of correct size
    set.seed(x)
    localSampleSize <- nrRepList[[x]]
    
    if(TRUE) {
      localSubset <- sample(
        colnames(countData),
        localSampleSize * 2
      )
      
      localDesign <- data.frame(
        sample_id = localSubset,
        condition = c(
          rep('a',floor  ( length(localSampleSize) )),
          rep('b',ceiling( length(localSampleSize) ))
        ),
        stringsAsFactors = FALSE
      )
      localDesign <- localDesign[sort.list(localDesign$condition),]
    }
    
    # Step 2: Subset to expressed features using edgeR::filterByExpr or DRIMSeq::dmFilter
    if(filter=="edgeR") {
      
      y <- edgeR::DGEList(counts = countData[,localDesign$sample_id])
      design <- model.matrix(~condition, data=localDesign)
      
      filter <- edgeR::filterByExpr(y,design=design,
                                    min.count = edgeR_filter_spec$min.count, 
                                    min.total.count = edgeR_filter_spec$min.total.count, 
                                    large.n = edgeR_filter_spec$large.n, 
                                    min.prop = edgeR_filter_spec$min.prop) 
      
      localCm <- y$counts[filter,]
      # Get only multi-isoform genes (after filtering)
      localTx <- metaData[which(
        metaData$isoform_id %in% rownames(localCm)),]
      
      tmp <- table(localTx$gene_id)
      tmp <- tmp[which( tmp >= 2)]
      
      localTx <- localTx[which(localTx$gene_id %in% names(tmp)),]
      localCm <- localCm[which(rownames(localCm) %in% localTx$isoform_id),]
    } else if(filter=="DRIMSeq") {
      
      localCm <-  as.data.frame(countData[,localDesign$sample_id])
      
      geneForEachTx <- metaData[match(rownames(localCm),metaData[,"isoform_id"]),"gene_id"]
      
      localCm$gene_id <- geneForEachTx
      localCm$feature_id <- row.names(localCm)
      
      d <- DRIMSeq::dmDSdata(counts = localCm, samples = localDesign)

      d_filter <- dmFilter(d,
                           min_samps_feature_expr=localSampleSize/2, 
                           min_feature_expr=10, 
                           min_samps_feature_prop=localSampleSize/2, 
                           min_feature_prop=0.1,
                           min_samps_gene_expr=localSampleSize, 
                           min_gene_expr=10)
      
      localCm <- countData[counts(d_filter)$feature_id,localDesign$sample_id]
      
      # Get only multi-isoform genes (after filtering)
      localTx <- metaData[which(
        metaData$isoform_id %in% rownames(localCm)),]
      
      tmp <- table(localTx$gene_id)
      tmp <- tmp[which( tmp >= 2)]
      
      localTx <- localTx[which(localTx$gene_id %in% names(tmp)),]
      localCm <- localCm[which(rownames(localCm) %in% localTx$isoform_id),]
    }
    
    # Step 3: Extract isoforms to modify and swap their expressions
    if(TRUE) {
      genesToModify <- sample(
        x = unique(localTx$gene_id),
        size = round(
          length(unique(localTx$gene_id)) * fracGenesAffected
        )
      )
      
      samplesToModify <- localDesign$sample_id[which(
        localDesign$condition == 'b')]
      
      cm_swapped <- localCm
      
      transcripts_toSwap_current <- c()
      transcripts_swapped_current <- c()
      transcripts_swapped_all <- c()
      transcripts_toSwap_all <- c()
      
      for (gene in genesToModify) {
        
        current <- localTx[which(localTx$gene_id==gene),]
        nSwap <- max(2,rbinom(1,nrow(current),1/3))
        
        transcripts_toSwap_current <- sample(
          x = current$isoform_id,
          size = nSwap)
        
        # swap order of txs completely by putting the fist one last
        transcripts_swapped_current <- c(transcripts_toSwap_current[-1],transcripts_toSwap_current[1])
        
        # Add to swapping queue
        transcripts_swapped_all <- c(transcripts_swapped_all,transcripts_swapped_current)
        transcripts_toSwap_all <- c(transcripts_toSwap_all,transcripts_toSwap_current)
      }
      
      # Perform the swapping in matrix
      cm_swapped[transcripts_toSwap_all,which(colnames(cm_swapped)%in%samplesToModify)] <- cm_swapped[transcripts_swapped_all,which(colnames(cm_swapped)%in%samplesToModify)]
      
      # Set swapping in localTx
      localTx$txSwapped <- vector(length=nrow(localTx))
      localTx[which(localTx$isoform_id %in% transcripts_swapped_all),"txSwapped"] <- TRUE 
      localTx$nrSamplesPerCondition <- localSampleSize
    }
    
    colnames(localTx) <- c('TXNAME','GENEID','txSwapped','nrSamplesPerCondition')
    
    # Combine data
    dataList <- list(
      data     = cm_swapped,
      design   = localDesign,
      metaInfo = localTx
    )
    
    return(dataList)
  })
}
