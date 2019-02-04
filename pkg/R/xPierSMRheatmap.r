#' Function to visualise SMR evidence using heatmap
#'
#' \code{xPierSMRheatmap} is supposed to visualise SMR evidence using heatmap. It returns an object of class "ggplot".
#'
#' @param data an input vector containing gene symbols
#' @param xTarget a "dTarget" or "sTarget" object with the componet 'list_pNode' related to 'eGene' predictors. Alternatively, it can be a data frame with columns ('Context','mode','probeID','Gene','ProbeChr','Probe_bp','topSNP','topSNP_chr','topSNP_bp','A1','A2','b_GWAS','b_eQTL','b_SMR','p_GWAS','p_eQTL','p_SMR','fdr_SMR')
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param zlim the minimum and maximum z values for which colors should be plotted
#' @return an object of class "ggplot" appended with 'mat' (the matrix colored by 'b_SMR') and 'df' (a data frame with columns 'priority','code','Context','mode','probeID','Gene','ProbeChr','Probe_bp','topSNP','topSNP_chr','topSNP_bp','A1','A2','b_GWAS','b_eQTL','b_SMR','p_GWAS','p_eQTL','p_SMR','fdr_SMR','direction_GWAS','direction_eQTL','direction_SMR'). Only those fdr_SMR<0.05 are visualised.
#' @note none
#' @export
#' @seealso \code{\link{xPierSMRheatmap}}
#' @include xPierSMRheatmap.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata/"
#' 
#' gp <- xPierSMRheatmap(data, dTarget)
#' }

xPierSMRheatmap <- function(data, xTarget, colormap='steelblue-white-darkred', zlim=c(-1,1))
{
    
    ## df_evidence
    if(class(xTarget) == "sTarget" | class(xTarget) == "dTarget"){
    
    	if(any(names(xTarget) %in% 'list_pNode')){
			list_pNode <- xTarget$list_pNode
			ind <- which(grepl('eGene_', names(list_pNode)))
			ls_df <- lapply(list_pNode[ind], function(x){
				x$evidence
			})
			df_evidence <- do.call(rbind, ls_df)
		
			## df_gene_priority: extract priority
			df_gene_priority <- data.frame(gene=xTarget$priority$name, priority=xTarget$priority$priority, stringsAsFactors=F)
			ind <- match(df_gene_priority$gene, data)
			if(sum(!is.na(ind))>0){
				df_gene_priority <- df_gene_priority[!is.na(ind),]
			}
		
		}else{
			return(NULL)
		}
		
    }else if(class(xTarget) == "data.frame"){
    	df_evidence <- xTarget[,c('Context','mode','probeID','Gene','ProbeChr','Probe_bp','topSNP','topSNP_chr','topSNP_bp','A1','A2','b_GWAS','b_eQTL','b_SMR','p_GWAS','p_eQTL','p_SMR','fdr_SMR')]
    	
		## df_gene_priority: always priority=1
		df_gene_priority <- data.frame(gene=data, priority=1, stringsAsFactors=F)
    	
    }else{
    	return(NULL)
    	
    }
    
    ######################################################################################
    
	## df_evidence_priority
	ind <- match(df_evidence$Gene, df_gene_priority$gene)
	df_evidence_priority <- df_evidence[!is.na(ind),]
	df_evidence_priority$priority <- df_gene_priority$priority[ind[!is.na(ind)]]

	######################################################################################
	
	Context <- code <- p_SMR <- priority <- fdr_SMR <- NULL
	
	## append 'direction_GWAS', 'direction_eQTL', 'direction_SMR' and 'code'
	df <- df_evidence_priority
	df$direction_GWAS <- ifelse(df$b_GWAS>0, 'RSK', 'PRT')
	df$direction_eQTL <- ifelse(df$b_eQTL>0, 'ASC', 'DES')
	df$direction_SMR <- ifelse(df$b_SMR>0, 'POS', 'NEG')
	df$code <- paste0(df$Gene,' (',df$topSNP,') [',df$A1,'<',df$A2,'; ',df$direction_GWAS,'; ',df$direction_eQTL,']')
	
	## keep the best p_SMR per Context and code
	df <- df %>% dplyr::arrange(Context,code,p_SMR)
	ind <- which(!duplicated(df[,c("Context","code")]))
	df <- df[ind,]	
	df <- subset(df, fdr_SMR < 0.05)
	df <- df %>% dplyr::arrange(-priority)
	
	## columns alwawys sorted by pre-defined ones
	## rows always sorted by priority
    default.Context <- c("CD14","LPS2","LPS24","IFN","Bcell","NK","Neutrophil","CD4","CD8","Blood","Monocyte","shared_CD14","shared_LPS2","shared_LPS24","shared_IFN")
	ind <- match(default.Context, unique(df$Context))
	columns <- default.Context[!is.na(ind)]
	mat <- as.matrix(xSparseMatrix(df[,c('code','Context','b_SMR')], rows=unique(df$code), columns=columns, verbose=F))
	mat[mat==0] <- NA
    gp <- xHeatmap(mat, reorder="none", colormap=colormap, ncolors=64, zlim=zlim, legend.title="Effect", barwidth=0.4, x.rotate=60, shape=19, size=2, x.text.size=7, y.text.size=6, legend.text.size=5, legend.title.size=7, na.color='transparent', barheight=max(3,min(5,nrow(mat))))
    
    ## append 'mat' and 'df'
	gp$mat <- mat
	gp$df <- df[,c('priority','code','Context','mode','probeID','Gene','ProbeChr','Probe_bp','topSNP','topSNP_chr','topSNP_bp','A1','A2','b_GWAS','b_eQTL','b_SMR','p_GWAS','p_eQTL','p_SMR','fdr_SMR','direction_GWAS','direction_eQTL','direction_SMR')]
	
	invisible(gp)
}

