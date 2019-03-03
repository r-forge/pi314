#' Function to prioritise genes based on seed eGenes identified through SMR integrating GWAS and eQTL summary data
#'
#' \code{xPierSMR} is supposed to prioritise genes based on seed eGenes identified through SMR integrating GWAS and eQTL summary data. To prioritise genes, it first conducts Summary-data-based Mendelian Randomisation (SMR) integrating GWAS and eQTL summary data to identify and score seed genes (that is, eGenes weighted by SMR reported p-values). It implements Random Walk with Restart (RWR) and calculates the affinity score of all nodes in the graph to the seeds. The priority score is the affinity score. Parallel computing is also supported for Linux-like or Windows operating systems. It returns an object of class "pNode". 
#'
#' @param data a data frame storing GWAS summary data with following required columns 'snp', 'effect' (the effect allele assessed), 'other' (other allele), 'b' (effect size for the allele assessed; log(odds ratio) for a case-control study), 'se' (standard error), 'p' (p-value), and the optional columnns 'freq' (frequency of the effect allele; not essential unless 'freq.check' is true) and 'n' (sample size; not required)
#' @param eqtl context-specific eQTL summary data. It can be one of "Bcell","Blood","CD14","CD4","CD8","IFN","LPS24","LPS2","Monocyte","Neutrophil","NK","shared_CD14","shared_IFN","shared_LPS24","shared_LPS2"
#' @param peqtl eQTL p-value threshold for selecting a probe (with the top associated eQTL passing a p-value threshold) for the SMR test. In other words, a probe with the top associated eQTL not passing this threshold will be removed for the test. By default, it is 5e-2
#' @param heidi logical to indicate whether the HEIDI test is enabled. By default it is true
#' @param bfile a character specifying where to find the LD reference data containing three files (.bed, .bim, and .fam). Required if heidi test is enabled (see above)
#' @param clear logical to indicate whether the temporary and log files are cleared up. By default, it sets to TRUE
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathways Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addtion to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways
#' @param STRING.only the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used
#' @param weighted logical to indicate whether edge weights should be considered. By default, it sets to false. If true, it only works for the network from the STRING database 
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network. If the user provides the "igraph" object with the "weight" edge attribute, RWR will assume to walk on the weighted network
#' @param seeds.inclusive logical to indicate whether non-network seed genes are included for prioritisation. If TRUE (by default), these genes will be added to the netowrk
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for Random Walk with Restart (RWR). The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param parallel logical to indicate whether parallel computation with multicores is used. By default, it sets to true, but not necessarily does so. Partly because parallel backends available will be system-specific (now only Linux or Mac OS). Also, it will depend on whether these two packages "foreach" and "doMC" have been installed
#' @param multicores an integer to specify how many cores will be registered as the multicore parallel backend to the 'foreach' package. If NULL, it will use a half of cores available in a user's computer. This option only works when parallel computation is enabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param ... additional graphic parameters used in xMEsmr. They are "mode", "window.cis", "window.trans", "heidi.peqtl", "heidi.ld", "heidi.num", "freq.check", "thread.num", "p.adjust.method" and "silent"
#' @return
#' an object of class "pNode", a list with following components:
#' \itemize{
#'  \item{\code{priority}: a matrix of nNode X 6 containing node priority information, where nNode is the number of nodes in the input graph, and the 5 columns are "name" (node names), "node" (1 for network genes, 0 for non-network seed genes), "seed" (1 for seeds, 0 for non-seeds), "weight" (weight values),  "priority" (the priority scores that are rescaled to the range [0,1]), "rank" (ranks of the priority scores), "description" (node description)}
#'  \item{\code{g}: an input "igraph" object}
#'  \item{\code{evidence}: a data frame storing evidence}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note The input graph will treat as an unweighted graph if there is no 'weight' edge attribute associated with
#' @export
#' @seealso \code{\link{xPierGenes}}
#' @include xPierSMR.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' data <- utils::read.delim(file="summary_gwas.RA.txt", header=T, row.names=NULL, stringsAsFactors=F)
#' pNode_smr <- xPierSMR(data, eqtl="Blood", network="STRING_high", restart=0.7, RData.location=RData.location)
#' write.table(pNode_smr$priority, file="Genes_priority.SMR.txt", sep="\t", row.names=FALSE)
#' }

xPierSMR <- function(data, eqtl=c("CD14","LPS2","LPS24","IFN","Bcell","NK","Neutrophil","CD4","CD8","Blood","Monocyte","shared_CD14","shared_LPS2","shared_LPS24","shared_IFN"), peqtl=5e-8, heidi=F, bfile=NULL, clear=T, network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD","KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], weighted=FALSE, network.customised=NULL, seeds.inclusive=TRUE, normalise=c("laplacian","row","column","none"), restart=0.7, normalise.affinity.matrix=c("none","quantile"), parallel=TRUE, multicores=NULL, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata", ...)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    eqtl <- match.arg(eqtl)
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    pNode <- NULL
    se <- b <- effect <- other <- n <- p <- NULL
	Gene <- p_SMR <- fdr_SMR <- NULL
	
	###########################	
	# summary GWAS
	###########################
	if(class(data)=='data.frame'){
		if(all(c('snp','effect','other','b','p','se') %in% colnames(data))){
			df <- data[,c('snp','effect','other','b','p','se')]
			
			if('freq' %in% colnames(data)){
				df$freq <- data$freq
			}else{
				df$freq <- NA
			}
			if('n' %in% colnames(data)){
				df$n <- data$n
			}else{
				df$n <- NA
			}
			
			############
			# random numbers derived from time stamps
			tempnum <- gsub('\\.*', '', as.character(as.numeric(Sys.time())))
			gwas.summary <- paste0(tempnum,'_gwas.txt')
			############
			output <- df[,c("snp","effect","other","freq","b","p","se")] %>% dplyr::filter(!is.na(se) & b!=0 & effect!='' & other!='') %>% dplyr::mutate(n=NA) %>% dplyr::filter(p<5e-8)
			utils::write.table(output, file=gwas.summary, row.names=F, col.names=T, quote=F, sep="\t")
				
		}else{
			warnings("The input data.frame does not contain required columns: c('snp','effect','other','b','p','se').\n")
			return(NULL)
		}
	}
		
	###########################	
	# summary eQTL
	###########################
	## beqtl.summary
	beqtl.summary <- file.path(RData.location, "Pi_eQTL_hg19", eqtl)
	vec <- paste0(beqtl.summary, c('.besd','.esi','.epi'))
	if(any(!file.exists(vec))){
		if(verbose){
			message(sprintf("The beqtl.summary '%s' not found (%s)!", beqtl.summary, as.character(Sys.time())), appendLF=T)
		}
			
		##############
		## create a new directory to hold the downloads
		my_dir <- file.path(getwd(), "Pi_eQTL_hg19")
		if(verbose){
			message(sprintf("\tcreate directory: '%s' (%s) ...", my_dir, as.character(Sys.time())), appendLF=T)
		}
		if(!file.exists(my_dir)){
			dir.create(my_dir)
		}
		if(verbose){
			message(sprintf("\tdownloading files (once and only once) into '%s' (%s) ...", my_dir, as.character(Sys.time())), appendLF=T)
		}
		vec_files <- paste0(eqtl, c('.besd','.esi','.epi'))
		source_files <- file.path("http://galahad.well.ox.ac.uk/bigdata","Pi_eQTL_hg19",vec_files)
		## download all files
		ls_tmp <- lapply(source_files, function(x){
			source <- x
			target <- file.path(my_dir, basename(x))
			utils::download.file(source, target, quiet=T)
		})
		##############
		beqtl.summary <- file.path(my_dir, eqtl)
		if(verbose){
			message(sprintf("The beqtl.summary '%s' used instead (%s)!", beqtl.summary, as.character(Sys.time())), appendLF=T)
		}
		
	}

	###########################	
	# bfile
	###########################
	if(heidi & is.null(bfile)){
		## bfile
		bfile <- file.path(RData.location, "Pi_eQTL_hg19", "Merged_EUR")
		vec <- paste0(bfile, c('.bed','.fam','.bim'))
		if(any(!file.exists(vec))){
			if(verbose){
				message(sprintf("The bfile '%s' not found (%s)!", bfile, as.character(Sys.time())), appendLF=T)
			}
			##############
			## create a new directory to hold the downloads
			my_dir <- file.path(getwd(), "Pi_eQTL_hg19")
			if(verbose){
				message(sprintf("\tcreate directory: '%s' (%s) ...", my_dir, as.character(Sys.time())), appendLF=T)
			}
			if(!file.exists(my_dir)){
				dir.create(my_dir)
			}
			if(verbose){
				message(sprintf("\tdownloading files (once and only once) into '%s' (%s) ...", my_dir, as.character(Sys.time())), appendLF=T)
			}
			vec_files <- paste0("Merged_EUR", c('.bed','.fam','.bim'))
			source_files <- file.path("http://galahad.well.ox.ac.uk/bigdata","Pi_eQTL_hg19",vec_files)
			## download all files
			ls_tmp <- lapply(source_files, function(x){
				source <- x
				target <- file.path(my_dir, basename(x))
				utils::download.file(source, target, quiet=T)
			})
			##############
			bfile <- file.path(my_dir, "Merged_EUR")
			if(verbose){
				message(sprintf("The bfile '%s' used instead (%s)!", bfile, as.character(Sys.time())), appendLF=T)
			}
		
		}
	}
	
	#df_output <- xMEsmr(gwas.summary, beqtl.summary, peqtl=peqtl, heidi=heidi, bfile=bfile, clear=clear, verbose=verbose, ...)
	if(class(suppressWarnings(try(df_output <- xMEsmr(gwas.summary, beqtl.summary, peqtl=peqtl, heidi=heidi, bfile=bfile, clear=clear, verbose=verbose, ...), T)))=="try-error"){
		df_output <- NULL
	}
	
	if(class(df_output)=='data.frame'){
		
		if(0){
			#################################
			# remove HLA genes and histone genes
			ind <- which(!grepl('^HLA-|^HIST', df_output$Gene))
			df_output <- df_output[ind,]
			#################################
		}
		
		df_evidence <- data.frame(Context=eqtl, df_output, stringsAsFactors=F)
		
		###########
		df_output <- df_output %>% dplyr::arrange(Gene,p_SMR)
		ind <- which(!duplicated(df_output[,c("Gene")]))
		df_output <- df_output[ind,]
		
		## pass SMR test: fdr_SMR<0.05
		df_output <- subset(df_output, fdr_SMR<0.05)
		## also pass HEIDI test: fdr_HEIDI>=0.05
		if(heidi){
			fdr_HEIDI <- NULL
			df_output <- subset(df_output, fdr_HEIDI>=0.05 | is.na(fdr_HEIDI))
		}
		
		## the seed gene (eGene) weighted by p_SMR
		data_subset <- df_output[,c("Gene","p_SMR")]
		data_subset$p_SMR <- -log10(data_subset$p_SMR)
		###########
		
    	pNode <- suppressMessages(xPierGenes(data=data_subset, network=network, STRING.only=STRING.only, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose, RData.location=RData.location))
    	
    	pNode$evidence <- df_evidence
    	
		if(clear){
			## remove gwas.summary (*_gwas.txt)
			unlink(gwas.summary)
		}
    }
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    }
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(pNode)
}
