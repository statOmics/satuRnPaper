suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(BANDITS))

run_BANDITS <- function(L,countData,tx2gene) {
  message("BANDITS")
  session_info <- sessionInfo()
  timing <- system.time({
    
    tx2gene <- tx2gene_init
    quantsf_counts <- countData_init
    removed <- "none"

    quantsf_counts <- quantsf_counts[,-which(colnames(quantsf_counts) %in% c("SRR2727161","SRR2727062","SRR2727475","SRR2727100"))] 
    # ECC failed for these samples
    
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

      	## get ECC files
      	dir <- "/home/compomics/salmon_Chen_new/"
      	samples <- read.table(file.path(dir, "SRRIDChen.txt"), header = FALSE)

	dir <- "/home/compomics/salmon_Chen_new"      
      	equiv_classes_files = file.path(dir, as.character(samples$V1), "aux_info", "eq_classes.txt")
      	equiv_classes_files <- equiv_classes_files[which(substring(equiv_classes_files,33,42)%in%colnames(quantsf_counts))]
	
	message(equiv_classes_files[!file.exists(equiv_classes_files)])

		
	eff_len <- readRDS("/home/compomics/scTime_Chen/eff_len.Rds")

	eff_len <- eff_len[which(names(eff_len) %in% tx2gene$TXNAME)]

	tx2gene <- tx2gene[,c(2,1)]

	transcripts_to_keep <- filter_transcripts(gene_to_transcript = tx2gene,
                                          transcript_counts = quantsf_counts,
                                          min_transcript_proportion = 0,
                                          min_transcript_counts = 0,
                                          min_gene_counts = 0)

	rm(list=setdiff(ls(),c("transcripts_to_keep","quantsf_counts","tx2gene","design","equiv_classes_files","eff_len","timing","session_info","removed")))
	
	print(Sys.time())
	input_data <- create_data(salmon_or_kallisto = "salmon",
                         gene_to_transcript = tx2gene,
                         salmon_path_to_eq_classes = equiv_classes_files,
                         eff_len = eff_len, 
                         n_cores = 1,
                         transcripts_to_keep = transcripts_to_keep)
	print(Sys.time())

	precision <- prior_precision(gene_to_transcript = tx2gene,
                            transcript_counts = quantsf_counts,
                            n_cores = 1)
	
	print(Sys.time())
	results <- BANDITS::test_DTU(BANDITS_data = input_data,
                   precision = precision$prior,
                   samples_design = design,
                   group_col_name = "group",
                   R = 10^4, 
                   burn_in = 2*10^3,
                   gene_to_transcript = tx2gene,
                   n_cores = 1)
	print(Sys.time())

	finalResult <- nrow(results@Transcript_results[which(results@Transcript_results$adj.p.values < 0.05),])

	mem <- gc()
    
    } else {
      ## make sure to have an output for result and info
        quantsf_counts <- matrix(data=NA, nrow=1,ncol=1)
        finalResult <- "More TXs asked than there were present"
	
	mem <- gc()
    }
  })

  output_BANDITS <- list(session_info = session_info,info = list(size = i,TXs = nrow(quantsf_counts)),timing = timing,removed = removed,memory = mem,result = finalResult)
  
  assign("result", output_BANDITS, envir=globalenv())
  rm(list = ls())
  gc()
}
