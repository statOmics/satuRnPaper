suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(locfdr))
suppressPackageStartupMessages(library(speedglm))

## get all the required functions first (analogous to load  library, but library doesn't exist yet)

# Helper functions for fitting the quasibinomial

StatModel<-function(type="fitError",params=list(),varPosterior=is.numeric(NA),dfPosterior=is.numeric(NA)) {
    out<-new("StatModel")
    out@type=type
    out@params=params
    out@varPosterior=varPosterior
    out@dfPosterior=dfPosterior
    return(out)
}

.StatModel <- setClass("StatModel",
                        slots = c(type = "character",
                        params = "list",
                        varPosterior = "numeric",
                        dfPosterior = "numeric"))

getModel <- function(object) {
  return(object@params)
}

calcDispersion <- function(model,type){

    if(type != "fitError") {
        df.r <- model$df.residual
        if(df.r > 0){
            #if(any(mod$weights==0))
		        #warning("observations with zero weight not used for calculating dispersion")
		        model$dispersion <- sum((model$weights*model$residuals^2)[model$weights > 0])/ df.r
        } else {
            model$dispersion <- NA
        }
    } else {
        model$dispersion <- NA
    }
    return(model)
}

getDispersion <- function(object){
    if(object@type!="fitError") {
        dispersion <- object@params$dispersion
    } else {
        dispersion <- NA
    }
    return(dispersion)
}

getDF <- function(object){
    if (!object@type=="fitError") {
      df <- object@params$df.residual
      #df <- object@params$df
    } else {
      df <- NA
    }
    return(df)
}

vcovUnscaled <- function(model,type) {
  
    if (!type=="fitError") {
      p1 <- 1L:model$rank
      p <- length(model$coef)
      out <- matrix(NA,p,p)
      out[p1,p1] <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
      colnames(out) <- rownames(out) <- names(model$coefficients)
      model$vcovUnsc <- out
    } else {
      model$vcovUnsc <- NA
    }
    return(model)
}

getOtherCount <- function(countData, tx2gene){
    # get tx2gene in better format
    geneForEachTx <- tx2gene$gene_id[match(rownames(countData),tx2gene$isoform_id)]
    geneForEachTx <- as.character(geneForEachTx)
    stopifnot(class(geneForEachTx) %in% c("character", "factor"))
    stopifnot(length(geneForEachTx) == nrow(countData))
    
    forCycle <- split(1:nrow(countData), as.character(geneForEachTx))
    all <- lapply(forCycle, function(i) {
        sct <- countData[i, , drop = FALSE]
        rs <- t(sapply(1:nrow(sct), function(r) colSums(sct[-r, , drop = FALSE])))
        rownames(rs) <- rownames(sct)
        rs
    })
        
    otherCount <- do.call(rbind, all)
    otherCount <- otherCount[rownames(countData), ]
    return(otherCount)
}


# Fit function for the quasibinomial

fit.qb <- function(countData, tx2gene, design, speed = FALSE){
  
    stopifnot(class(countData) %in% c("matrix", "data.frame"))
    countData <- as.matrix(countData)

    ## get the "other" counts
    otherCount <- getOtherCount(countData, tx2gene)
    stopifnot(all(rownames(countData) %in% rownames(otherCount)))
    
    ## The actual fit function
    fitQuasiLogistic <- function(countData, otherCount, design, speed = speed){
        
        countsAll <- cbind(countData, otherCount)
        drop <- rowSums(countsAll) == 0 ## gene count is zero in this cell
        countsAll <- countsAll + 1
        countsAll[drop,] <- NA

        if (speed == FALSE) {
          model <- try(glm(countsAll ~ -1+design, family="quasibinomial"))
        } else {
          model <- try(speedglm(countsAll ~ -1+design, family=quasibinomial()))
        }
        #model <- try(speedglm(countsAll ~ -1+design, family=quasibinomial()))
        #model <- try(fastglm(design, counts, family=quasibinomial(), method=3, maxit = 25))
        if (class(model)[1]=="try-error"){ 
            # catch error
            model <- list()
            type <- "fitError"
            class(model) <- "list"
        } else {
            type <- "glm"
            class(model) <- "list"
        }
        
        if (speed == FALSE){
          model <- calcDispersion(model,type) ## calculate disp slot
          model <- vcovUnscaled(model,type) ## calculate vcov slot
          model <- model[c("coefficients","df.residual", "dispersion", "vcovUnsc")]
        } else {
          if (type != "fitError") {
            model$vcovUnsc <- solve(model$XTX)
          } else {
            model$vcovUnsc <- NA
          }
          
          model <- model[c("coefficients","df", "dispersion", "vcovUnsc")]
          names(model) <- c("coefficients","df.residual", "dispersion", "vcovUnsc")
        }

        .out  <- .StatModel(type = type, 
                      params = model,
                      varPosterior = as.numeric(NA),
                      dfPosterior = as.numeric(NA))
        return(.out)
    }
    
    models <- pblapply(seq_len(nrow(countData)), function(i) fitQuasiLogistic(countData=countData[i,], otherCount=otherCount[i,], design = design, speed = speed))
    # retain transcript names 
    names(models) <- rownames(countData)

    # Squeeze a set of sample variances together by computing empirical Bayes posterior means
    hlp <- limma::squeezeVar(var = sapply(models, getDispersion),
                           df = sapply(models, getDF))
  
    # put variance and degrees of freedom in appropriate slots
    for (i in 1:length(models)) {
        mydf <- hlp$df.prior + getDF(models[[i]])
        models[[i]]@varPosterior <- as.numeric(hlp$var.post[i])
        models[[i]]@dfPosterior <- as.numeric(mydf)
    }

    # return object of class StatModel
    return(models)
}

# Helper functions for the topTable (test) function

getCoef <- function(object,L) {
  
  if (!is.null(object@params$coefficients)){
    coef <- object@params$coefficients
  } else {
    coef <- rep(NA, times=nrow(L))
  }
  return(coef)
}

getEstimates <- function(object,L){
  coefs <- getCoef(object,L)
  return(L%*%coefs)
}

varContrast <- function(object,L){
  if (object@type!="fitError"){
      if (nrow(object@params$vcovUnsc) == length(L)){
          vcovTmp <- object@params$vcovUnsc*object@varPosterior
          return(diag(t(L)%*%vcovTmp%*%L))
      }
  }
  return(NA)
}

getDfPosterior <- function(object){
  return(object@dfPosterior)
}

p.adjust_empirical_hlp <- function(zz){
  
  N <- length(zz)
  b <- 4.3 * exp(-0.26*log(N,10))
  med <- median(zz)
  sc <- diff(quantile(zz)[c(2,4)])/(2*qnorm(.75))
  mlests <- locfdr:::locmle(zz, xlim=c(med, b*sc)) ## default d=0, s=1 
  
  nulltype = 1
  lo <- min(zz)
  up <- max(zz)
  bre = 120
  breaks <- seq(lo, up, length = bre)
  zzz <- pmax(pmin(zz, up), lo)
  x <- (breaks[-1] + breaks[ - length(breaks)])/2 ## midpoints hist
  sw <- 0
  X <- cbind(1, poly(x, df = 7)) ## design
  zh <- hist(zzz, breaks = breaks, plot = F)
  y <- zh$counts
  f <- glm(y ~ poly(x, df = 7), poisson)$fit ## fitted values
  
  Cov.in = list(x=x, X=X, f=f, sw=sw)
  ml.out = locfdr:::locmle(zz, xlim = c(mlests[1], b * mlests[2]),
              d=mlests[1], s=mlests[2], Cov.in=Cov.in) ## try with without covin
  mlests = ml.out$mle
  
  return(mlests)
  
}

p.adjust_empirical <- function(pvalues,tvalues,plot=FALSE){
   zvalues <- qnorm(pvalues/2)*sign(tvalues)
   zvalues_mid <- zvalues[abs(zvalues) < 10] ## to avoid numerical issues
   zvalues_mid <- zvalues_mid[!is.na(zvalues_mid)]
   
   #locfdr_zval_null <- locfdr(zvalues_mid, nulltype = 1, plot=T)
   mlests <- p.adjust_empirical_hlp(zvalues_mid)

   zval_empirical <- (zvalues-mlests[1])/mlests[2]
   ## rescale with pi0
   pval_empirical <- 2*pnorm(-abs(zval_empirical), mean=0, sd=1)
   
   if (plot){
     
      zval_empirical <- zval_empirical[!is.na(zval_empirical)]
      zval_empirical <- zval_empirical[!is.infinite(zval_empirical)]
      lo <- min(zval_empirical)
      up <- max(zval_empirical)
      
      lo <- min(lo, -1*up)
      up <- max(up, -1*lo)
      
      bre = 120
      breaks <- seq(lo, up, length = bre)
      zzz <- pmax(pmin(zval_empirical, up), lo)
      zh <- hist(zzz, breaks = breaks, plot = F)
	    yall <- zh$counts
	    K <- length(yall)
     
      hist(zzz, breaks = breaks, xlab = " ", main = " ", freq = F)
		  xfit <- seq(min(zzz), max(zzz), length = 4000) 
      yfit <- dnorm(xfit, mean = 0, sd = 1)
      lines(xfit, yfit, col = "darkgreen", lwd = 2)
   }
   
   FDR <- p.adjust(pval_empirical, method = "BH")
   newList <- list("pval" = pval_empirical, "FDR" = FDR)
   ## also return pval_empirical
   return(newList)
}

# The topTable (test) function

topTable<-function(models,contrast, plot) {
    estimates <- sapply(models,getEstimates,L=contrast)
    se <- sqrt(sapply(models,varContrast,L=contrast)) # uses the squeezed variances
    df <- sapply(models,getDfPosterior)
    t <- estimates/se
    pval <- pt(-abs(t),df)*2
    regular_FDR <- p.adjust(pval, method = "BH") # regular FDR correction
    empirical <- p.adjust_empirical(pvalues=pval,tvalues=t,plot=plot) # empirical FDR correction
    empirical_pval <- empirical$pval
    empirical_FDR <-  empirical$FDR
    result <- data.frame(estimates,se,df,t,pval,regular_FDR,empirical_pval,empirical_FDR)
    #result <- result[order(result$empirical_FDR),]
    return(result)
}


## run function


run_quasiBinomial <- function(L,countData,tx2gene) {
  message("quasiBinomial")
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
        colnames(design) <- c("sampleId", "group")
        design$sampleId <- colnames(quantsf_counts)
        design$group <- sample(rep(c("A","B"),each=L$groupSizes[i]))
	
	## check if there are cells without counts due to filtering/subsampling
	## if so (unlikely); remove them, but store the fact they are removed as it may impact speed
	if(length(which(colSums(quantsf_counts) == 0)) > 0){
        	removed <- length(which(colSums(quantsf_counts) == 0))
        	design <- sampleData[-which(colSums(quantsf_counts) == 0),]
        	quantsf_counts <- quantsf_counts[,-which(colSums(quantsf_counts) == 0)]
      	}

	colnames(tx2gene)[1:2] <- c("isoform_id", "gene_id")
        group <- factor(design$group)
        design <- model.matrix(~0+group)
        
        models_glm <- fit.qb(countData = quantsf_counts, tx2gene = tx2gene, design = design, speed = FALSE)
        
        L <- matrix(0, ncol = 1, nrow = ncol(design))
        rownames(L) <- colnames(design)
        colnames(L) <- c("AvsB")
        L[1:2,] <- c(-1,1)
        
        res_glm <- apply(L, 2, function(i) topTable(models = models_glm, contrast = i, plot=FALSE))
        
        output <- data.frame(rownames(quantsf_counts))
        rownames(output) <- rownames(quantsf_counts)
        colnames(output) <- "TXNAME"
        
        output$GENEID <- tx2gene[match(output$TXNAME,tx2gene$TXNAME),"GENEID"]
        output$tvalue <- res_glm[[1]]$t
        output$p_value_raw <- res_glm[[1]]$pval
        output$FDR_raw <- res_glm[[1]]$regular_FDR
        output$p_value <- res_glm[[1]]$empirical_pval
        output$FDR <- res_glm[[1]]$empirical_FDR
    
        finalResult <- nrow(output[which(output$FDR < 0.05),])
	
	mem <- gc()
    
    } else {
      ## make sure to have an output for result and info
        quantsf_counts <- matrix(data=NA, nrow=1,ncol=1)
        finalResult <- "More TXs asked than there were present"
	
	mem <- gc()
    }
  })
  
  output_quasiBinomial <- list(session_info = session_info,
                                 info = list(size = i, 
                                             TXs = nrow(quantsf_counts)),
                                 timing = timing,
				 removed = removed,
				 memory = mem,
                                 result = finalResult)
  
  assign("result", output_quasiBinomial, envir=globalenv())
  rm(list = ls())
  gc()
}
