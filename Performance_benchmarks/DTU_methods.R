# Load all 7 DTU methods

satuRn_DTU <- function(countData, tx2gene, sampleData, quiet = FALSE){
    
    colnames(tx2gene)[1:2] <- c("isoform_id", "gene_id")
    group <- factor(sampleData$condition)
    design <- model.matrix(~0+group)
    
    sumExp <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=countData), colData = sampleData, rowData = tx2gene)
    
    metadata(sumExp)$formula <- ~0+group
    
    sumExp <- satuRn::fitQB(sumExp,
                            parallel = TRUE,
                            BPPARAM = BiocParallel::bpparam(),
                            verbose = TRUE
    )
    
    L <- matrix(0, ncol = 1, nrow = ncol(design))
    rownames(L) <- colnames(design)
    colnames(L) <- c("AvsB")
    L[1:2,] <- c(-1,1)
    
    sumExp <- satuRn::testDTU(object = sumExp, contrasts = L, plot=F)
    
    output <- data.frame(rownames(countData))
    rownames(output) <- rownames(countData)
    colnames(output) <- "TXNAME"
    
    output$GENEID <- tx2gene[match(output$TXNAME,tx2gene$TXNAME),"GENEID"]
    
    result <- rowData(sumExp)[["fitQBResult_AvsB"]]
    
    output$tvalue <- result$t
    output$p_value_raw <- result$pval
    output$FDR_raw <- result$regular_FDR
    output$p_value <- result$empirical_pval
    output$FDR <- result$empirical_FDR
    
    localRes <- output
    
    return(localRes)
}

## DoubleExpSeq
DoubleExpSeq_DTU <- function(countData, tx2gene, sampleData,quiet=FALSE) {
  
    # get tx2gene in better format
    geneForEachTx <- tx2gene$GENEID[match(rownames(countData),tx2gene$TXNAME)]
    groupID <- as.character(geneForEachTx)
    
    ## get the "total" counts
    stopifnot(class(countData)[1] %in% c("matrix", "data.frame"))
    countData <- as.matrix(countData)
    stopifnot(class(groupID) %in% c("character", "factor"))
    stopifnot(length(groupID) == nrow(countData))
    
    if( ! quiet) {print("Calculating offsets")}
    forCycle <- split(1:nrow(countData), as.character(groupID))
    all <- lapply(forCycle, function(i) {
      sct <- countData[i, , drop = FALSE]
      rs <- t(sapply(1:nrow(sct), function(r) colSums(sct[, , drop = FALSE])))
      rownames(rs) <- rownames(sct)
      rs
    })
    all <- do.call(rbind, all)
    stopifnot(all(rownames(countData) %in% rownames(all)))
    all <- all[rownames(countData), ]
    totalCount <- all
    
    condition <- as.factor(sampleData$condition)
    
    mDouble <- suppressWarnings(DBGLM1(y = countData,
                                       m = totalCount,
                                       groups = condition))
    
    pvalDouble <- mDouble$All[,"pVal"]
    dispDouble <- mDouble$All[,"Model.Disp"]
    
    localRes <- data.frame(rownames(countData))
    rownames(localRes) <- rownames(countData)
    colnames(localRes) <- "TXNAME"
    
    localRes$GENEID <- tx2gene[match(localRes$TXNAME,tx2gene$TXNAME),"GENEID"]
    localRes$p_value <- pvalDouble ## here I am only returning p-values, could also look into storing logFCs and FDRs
    localRes$disp <- dispDouble
    
    return(localRes)
}

### Limma diffsplice
limma_diffsplice_DTU <- function(countData, tx2gene, sampleData) {
    # Make design
    design = model.matrix(~condition, data=sampleData)
    
    geneForEachTx <- tx2gene$GENEID[match(rownames(countData),tx2gene$TXNAME)]
    y <- DGEList(counts = countData, group=as.data.frame(sampleData)$condition, genes = geneForEachTx)
    
    y <- calcNormFactors(y)
    
    v <- voom(y, design, plot=FALSE)
    fit <- lmFit(v, design)
    
    ex <- limma::diffSplice(fit, geneid = "genes")
    
    limmaRes <- topSplice(ex, coef=2, test="t", number = Inf)
    
    ### Massage result
    localRes <- limmaRes
    localRes$TXNAME <- rownames(localRes)
    rownames(localRes) <- NULL
    
    colnames(localRes)[which(colnames(localRes) == 'genes')] <- 'GENEID'
    
    localRes <- localRes[,c(
      which( colnames(localRes) == 'TXNAME'),
      which( colnames(localRes) == 'GENEID'),
      which( ! colnames(localRes) %in% c('TXNAME','GENEID'))
    )]
    
    colnames(localRes)[which(colnames(localRes) == 'P.Value')] <- 'p_value'
    
    return(localRes)
}

## edgeR_diffsplice
edgeR_diffsplice_DTU <- function(countData, tx2gene, sampleData) {
    
    y <- DGEList(counts = countData, group=as.data.frame(sampleData)$condition)
    y <- calcNormFactors(y)
    
    design <- model.matrix(~condition, data=sampleData)
    
    y <- estimateDisp(y, design = design)
    fit <- glmFit(y, design = design) ## glmQLfit --> very poor results
    
    tx2gene <- tx2gene[match(rownames(fit$counts), tx2gene$TXNAME),]
    
    dtuEdgeR <- edgeR::diffSpliceDGE(
      glmfit = fit,
      exonid = tx2gene$TXNAME,
      geneid = tx2gene$GENEID,
      coef = "conditionb",
      verbose = FALSE
    )
    
    localRes <- topSpliceDGE(
      dtuEdgeR,
      test = 'exon',
      number = Inf
    )
    
    rownames(localRes) <- NULL
    colnames(localRes)[which(colnames(localRes) == 'GeneID')] <- 'GENEID'
    colnames(localRes)[which(colnames(localRes) == 'ExonID')] <- 'TXNAME'
    
    localRes <- localRes[,c(
      which( colnames(localRes) == 'TXNAME'),
      which( colnames(localRes) == 'GENEID'),
      which( ! colnames(localRes) %in% c('TXNAME','GENEID'))
    )]
    
    colnames(localRes)[which(colnames(localRes) == 'P.Value')] <- 'p_value'
    
    return(localRes)
}

### DEXSeq
DEXSeq_DTU <- function(countData, tx2gene, sampleData) {
    # get tx2gene in better format
    geneForEachTx <- tx2gene$GENEID[match(rownames(countData),tx2gene$TXNAME)]
    
    countData <- round(countData)
    
    dxd <- DEXSeqDataSet(countData = countData,
                         sampleData = sampleData,
                         design = ~ sample + exon + condition:exon,
                         featureID = rownames(countData),
                         groupID = as.character(geneForEachTx))
    
    dxd <- estimateSizeFactors(dxd)
    ## note that estimating dispersions takes some time!
    print("estimating the dispersions")
    dxd <- estimateDispersions(dxd)
    
    print("testing for DEU")
    dxd <- testForDEU(dxd, reducedModel=~ sample + exon)
    
    print("Getting the results")
    dxr <- DEXSeqResults(dxd)
    
    ### Massage result
    localRes <- as.data.frame(dxr)
    localRes <- localRes[,c('featureID','groupID','exonBaseMean','dispersion','stat','pvalue','padj')]
    colnames(localRes)[1:2] <- c('TXNAME','GENEID')
    rownames(localRes) <- NULL
    
    colnames(localRes)[which(colnames(localRes) == 'pvalue')] <- 'p_value'
    colnames(localRes)[which(colnames(localRes) == 'padj')] <- 'FDR'
    
    return(localRes)
}

### DRIMSeq
DRIMSeq_DTU <- function(countData, tx2gene, sampleData) {
    ## To make DRIMSeq analysis reproducible
    set.seed(123)
    
    # get tx2gene in better format
    countData <- round(countData)
    countData <- as.data.frame(countData)
    
    geneForEachTx <- tx2gene$GENEID[match(rownames(countData),tx2gene$TXNAME)]
    
    countData$gene_id <- geneForEachTx
    countData$feature_id <- row.names(countData)
    rownames(countData) <- c()
    
    d <- dmDSdata(counts = countData, samples = sampleData)
    
    design <- model.matrix(~ condition, data = DRIMSeq::samples(d))
    
    print("dmPrecision")
    d <- dmPrecision(d, design = design)
    
    print("dmFit")
    d <- dmFit(d, design = design, verbose = 1)
    
    print("dmTest")
    d <- dmTest(d, coef = 2, verbose = 1)
    
    ### Massage result
    localRes <- DRIMSeq::results(d, level = "feature")[, c(
      'feature_id',
      'gene_id',
      'lr',
      'df',
      'pvalue',
      'adj_pvalue'
    )]
    colnames(localRes)[1:2] <- c('TXNAME','GENEID')
    
    colnames(localRes)[which(colnames(localRes) == 'pvalue')] <- 'p_value'
    colnames(localRes)[which(colnames(localRes) == 'adj_pvalue')] <- 'FDR'
    
    return(localRes)
}

## NBSplice
NBSplice_DTU <- function(countData, tx2gene, sampleData,quiet=FALSE) {
    
    sampleData$condition <- as.factor(sampleData$condition)
    geneIso <- as.data.frame(tx2gene[,c(2,1)])
    rownames(geneIso) <- geneIso[,2]
    colnames(geneIso) <- c("gene_id", "isoform_id")
    
    colnames(sampleData) <- c("id", "condition")
    rownames(sampleData) <- sampleData$id
    sampleData$condition <- as.factor(sampleData$condition)
    
    colName<-"condition"
    
    countData <- round(countData)
    
    myIsoDataSet <- IsoDataSet(countData, sampleData, colName, geneIso)
    
    myDSResults <- NBTest(myIsoDataSet, colName, test="F")
    

    localRes <- myDSResults@results
    localRes <- localRes[,c(1,2,7)]
    colnames(localRes) <- c("TXNAME", "GENEID", "p_value")
      
    return(localRes)
}
