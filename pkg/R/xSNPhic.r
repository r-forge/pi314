#' Function to extract promoter capture HiC-gene pairs given a list of SNPs
#'
#' \code{xSNPhic} is supposed to extract HiC-gene pairs given a list of SNPs.
#'
#' @param data NULL or a input vector containing SNPs. If NULL, all SNPs will be considered. If a input vector containing SNPs, SNPs should be provided as dbSNP ID (ie starting with rs) or in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is number; for example, 'chr16:28525386'. Alternatively, it can be other formats/entities (see the next parameter 'entity')
#' @param entity the data entity. By default, it is "SNP". For general use, it can also be one of "chr:start-end", "data.frame", "bed" or "GRanges"
#' @param include.HiC genes linked to input SNPs are also included. Genes linked to input SNPs are based on Promoter Capture HiC (PCHic), including 17 primary blood cell types: 'Monocytes','Macrophages_M0','Macrophages_M1','Macrophages_M2','Neutrophils','Megakaryocytes','Endothelial_precursors','Erythroblasts','Fetal_thymus','Naive_CD4_T_cells','Total_CD4_T_cells','Activated_total_CD4_T_cells','Nonactivated_total_CD4_T_cells','Naive_CD8_T_cells','Total_CD8_T_cells','Naive_B_cells','Total_B_cells'.
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' If input data is NULL, a data frame with following columns:
#' \itemize{
#'  \item{\code{from}: baited genomic regions (baits)}
#'  \item{\code{to}: preyed (other end) genomic regions of interactions (preys)}
#'  \item{\code{score}: CHiCAGO scores quantifying the strength of physical interactions between harbors and partners}
#' }
#' If input data is not NULL, a list with two components: "df" and "ig".
#' "df" is a data frame with following columns:
#' \itemize{
#'  \item{\code{from}: 'from/bait' genomic regions}
#'  \item{\code{to}: 'to/prey' genomic regions}
#'  \item{\code{score}: CHiCAGO scores quantifying the strength of physical interactions between baits and preys}
#'  \item{\code{from_genes}: genes associated with 'from/bait' genomic regions}
#'  \item{\code{to_genes}: genes associated with 'to/prey' genomic regions}
#'  \item{\code{SNP}: input SNPs (in query)}
#'  \item{\code{SNP_end}: specify which end SNPs in query fall into (either 'bait/from' or 'prey/to')}
#'  \item{\code{SNP_harbor}: genomic regions harbors the SNPs in query}
#'  \item{\code{Context}: the context in which PCHiC data was generated}
#' }
#' "ig" is an object of both classes "igraph" and "PCHiC", a direct graph with nodes for genomic regions and edges for CHiCAGO scores between them. Also added node attribute is 1) 'target' storing genes assocated and 2) 'SNP' for input SNPs (if the node harboring input SNPs). If several cell types are queried, "ig" is actually a list of "igraph"/"PCHiC" objects.
#' @note None
#' @export
#' @seealso \code{\link{xRDataLoader}}
#' @include xSNPhic.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' # a) provide the SNPs with the significance info
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' #data.file <- "http://galahad.well.ox.ac.uk/bigdata/AS.txt"
#' #AS <- read.delim(data.file, header=TRUE, stringsAsFactors=FALSE)
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' data <- names(ImmunoBase$AS$variants)
#'
#' \dontrun{
#' # b) extract HiC-gene pairs given a list of AS SNPs
#' PCHiC <- xSNPhic(data, include.HiC="Monocytes", GR.SNP="dbSNP_GWAS", RData.location=RData.location)
#' head(PCHiC$df)
#' 
#' # c) visualise the interaction (a directed graph: bait->prey)
#' g <- PCHiC$ig
#' ## a node with SNPs colored in 'skyblue' and the one without SNPs in 'pink'
#' ## the width in an edge is proportional to the interaction strength
#' xPCHiCplot(g, vertex.label.cex=0.5)
#' xPCHiCplot(g, glayout=layout_in_circle, vertex.label.cex=0.5)
#' }

xSNPhic <- function(data=NULL, entity=c("SNP","chr:start-end","data.frame","bed","GRanges"), include.HiC=c('Monocytes','Macrophages_M0','Macrophages_M1','Macrophages_M2','Neutrophils','Megakaryocytes','Endothelial_precursors','Erythroblasts','Fetal_thymus','Naive_CD4_T_cells','Total_CD4_T_cells','Activated_total_CD4_T_cells','Nonactivated_total_CD4_T_cells','Naive_CD8_T_cells','Total_CD8_T_cells','Naive_B_cells','Total_B_cells'), GR.SNP=c("dbSNP_GWAS","dbSNP_Common"), verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
	
	entity <- match.arg(entity)
	
    ######################################################
    # Link to targets based on HiC
    ######################################################
    
    default.include.HiC <- c('Monocytes','Macrophages_M0','Macrophages_M1','Macrophages_M2','Neutrophils','Megakaryocytes','Endothelial_precursors','Erythroblasts','Fetal_thymus','Naive_CD4_T_cells','Total_CD4_T_cells','Activated_total_CD4_T_cells','Nonactivated_total_CD4_T_cells','Naive_CD8_T_cells','Total_CD8_T_cells','Naive_B_cells','Total_B_cells')
	ind <- match(default.include.HiC, include.HiC)
	include.HiC <- default.include.HiC[!is.na(ind)]
    
    if(!is.null(data)){
    	if(entity=="SNP"){
    		data_gr <- xSNPlocations(data, GR.SNP=GR.SNP, verbose=verbose, RData.location=RData.location)
    	}else{
    		data_gr <- xGR(data, format=entity, verbose=verbose, RData.location=RData.location)
    	}
    	
    	if(is.null(data_gr)){
    		return(NULL)
    	}
    }
    
    df_returned <- NULL
    if(length(include.HiC) > 0){
    
		res_list <- lapply(include.HiC, function(x){

			if(verbose){
				now <- Sys.time()
				message(sprintf("Processing %s ...", x), appendLF=TRUE)
			}
			
			RData.customised <- paste('ig.PCHiC.', x, sep='')
			ig <- xRDataLoader(RData.customised=RData.customised, RData.location=RData.location, verbose=verbose)
			
			## Convert from igraph into data.frame
  			df_nodes <- igraph::get.data.frame(ig, what="vertices")
			df_edges <- igraph::get.data.frame(ig, what="edges")
			
			if(!is.null(data)){
			
				nodes_gr <- xGR(data=df_nodes[,1], format="chr:start-end", verbose=verbose, RData.location=RData.location)
				
				maxgap <- 0
				minoverlap <- 1L # 1b overlaps
				subject <- nodes_gr
				query <- data_gr
				q2r <- as.matrix(suppressWarnings(GenomicRanges::findOverlaps(query=query, subject=subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=TRUE)))
				
				res_df <- data.frame(SNP=names(data_gr)[q2r[,1]], nodes=names(nodes_gr)[q2r[,2]], stringsAsFactors=FALSE)
				
				ind_from <- match(res_df$nodes, df_edges[,'from'])
				nodes_from <- res_df[!is.na(ind_from),]
				partner_from <- df_edges[ind_from[!is.na(ind_from)], c('to','score')]
				partner_from_gene <- df_nodes[match(partner_from[,1],df_nodes[,1]), 2]
				nodes_from_gene <- df_nodes[match(nodes_from[,2],df_nodes[,1]), 2]
				df_from <- cbind(nodes_from, partner_from, nodes_from_gene, partner_from_gene, stringsAsFactors=FALSE)
				colnames(df_from) <- c('SNP', 'harbor', 'partner', 'score', 'harbor_genes', 'partner_genes')
				df_from$harbor_end <- rep('bait/from',nrow(df_from))
				
				ind_to <- match(res_df$nodes, df_edges[,'to'])
				nodes_to <- res_df[!is.na(ind_to),]
				partner_to <- df_edges[ind_to[!is.na(ind_to)], c('from','score')]
				partner_to_gene <- df_nodes[match(partner_to[,1],df_nodes[,1]), 2]
				nodes_to_gene <- df_nodes[match(nodes_to[,2],df_nodes[,1]), 2]
				df_to <- cbind(nodes_to, partner_to, nodes_to_gene, partner_to_gene, stringsAsFactors=FALSE)
				colnames(df_to) <- c('SNP', 'harbor', 'partner', 'score', 'harbor_genes', 'partner_genes')
				df_to$harbor_end <- rep('prey/to',nrow(df_to))
				
				df <- rbind(df_from, df_to)
				df$Context <- rep(x, nrow(df))
				
				####################
				# reverse conversion
				####################
				y <- df
				ind <- which(y$harbor_end=="bait/from")
				df_bait <- y[ind,c('harbor','partner','score','harbor_genes','partner_genes','SNP','harbor_end','harbor','Context')]
				colnames(df_bait) <- c('from','to','score','from_genes','to_genes','SNP','SNP_end','SNP_harbor','Context')
				ind <- which(y$harbor_end=="prey/to")
				df_prey <- y[ind,c('partner','harbor','score','partner_genes','harbor_genes','SNP','harbor_end','harbor','Context')]
				colnames(df_prey) <- c('from','to','score','from_genes','to_genes','SNP','SNP_end','SNP_harbor','Context')
				df <- rbind(df_bait, df_prey)
				rownames(df) <- NULL
				
			}else{
				df <- df_edges
				colnames(df) <- c('from', 'to', 'score')
				df$Context <- rep(x, nrow(df))
			}
			
			return(df)
		})
		## get data frame:
		### from to score Context
		### from to score from_genes to_genes SNP SNP_end SNP_harbor Context
		df_returned <- do.call(rbind, res_list)
	
	}
	
	if(!is.null(data)){
		if(verbose){
			now <- Sys.time()
			message(sprintf("Amongst %d SNPs, %d SNPs are falling into %d interacton regions", length(unique(data)), length(unique(df_returned$SNP)), length(unique(df_returned$SNP_harbor))), appendLF=TRUE)
		}
		
		####################################
		## also output igraph
		context_ls <- split(x=df_returned, f=df_returned$Context)
		ls_ig <- lapply(1:length(context_ls), function(i){
			df <- context_ls[[i]]
			
			if(1){
				#################################
				res_ls <- lapply(split(x=df$SNP, f=df$SNP_harbor),function(x){
					paste(x, collapse=';')
				})
				vec_harbor_SNPs <- unlist(res_ls)
				#################################
				
				## edges
				relations <- df[,1:3]
				relations <- relations[!duplicated(relations),]
				
				## nodes
				nodes <- base::as.data.frame(rbind(as.matrix(df[,c("from","from_genes")]), as.matrix(df[,c("to","to_genes")])), stringsAsFactors=FALSE)
				nodes <- nodes[!duplicated(nodes),]
				colnames(nodes) <- c("name","target")
				## apend list of SNP per harbor
				nodes$SNP <- NA
				ind <- match(nodes$name,names(vec_harbor_SNPs))
				nodes$SNP[!is.na(ind)] <- vec_harbor_SNPs[ind[!is.na(ind)]]

				ig <- graph.data.frame(d=relations, directed=TRUE, vertices=nodes)
				class(ig) <- c("PCHiC","igraph")
				return(ig)	
			}
			
		})
		names(ls_ig) <- names(context_ls)
		####################################
		
		if(length(ls_ig)==1){
			ls_ig <- ls_ig[[1]]
		}
		
		output <- list(df=df_returned, ig=ls_ig)
		
	}else{
		output <- df_returned
	}
    
    invisible(output)
}
