#' Function to extract priority matrix from a list of pNode objects
#'
#' \code{xPierMatrix} is supposed to extract priority matrix from a list of pNode objects. Also supported is the aggregation of priority matrix (similar to the meta-analysis) generating the priority results; we view this functionality as the discovery mode of the prioritisation.
#'
#' @param list_pNode a list of "pNode" objects or a "pNode" object
#' @param displayBy which priority will be extracted. It can be "score" for priority score (by default), "rank" for priority rank, "weight" for seed weight, "pvalue" for priority p-value
#' @param combineBy how to resolve nodes/targets from a list of "pNode" objects. It can be "intersect" for intersecting nodes (by default), "union" for unionising nodes
#' @param aggregateBy the aggregate method used. It can be either "none" for no aggregation, or "orderStatistic" for the method based on the order statistics of p-values, "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'fishers' and 'logistic' are preferred here
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' If aggregateBy is 'none' (by default), a data frame containing priority matrix, with each column/predictor for either priority score, or priorty rank or priority p-value.
#' If aggregateBy is not 'none', an object of the class "dTarget", a list with following components:
#' \itemize{
#'  \item{\code{priority}: a data frame of nGene X 6 containing gene priority (aggregated) information, where nGene is the number of genes, and the 6 columns are "name" (gene names), "rank" (ranks of the priority scores), "pvalue" (the aggregated p-value, converted from empirical cumulative distribution of the probability of being GSP), "fdr" (fdr adjusted from the aggregated p-value), "priority" (-log10(pvalue) but rescaled into the 5-star ratings), "description" (gene description) and seed info including "Overall" for the number of different types of seeds, followed by details on individual type of seeds (that is, "OMIM", "Phenotype", "Function", "nearbyGenes", "eQTL", "HiC")}
#'  \item{\code{predictor}: a data frame containing predictor matrix, with each column/predictor for either priority score, or priorty rank or priority p-value}
#'  \item{\code{metag}: an "igraph" object}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xPierSNPsAdv}}
#' @include xPierMatrix.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' df_score <- xPierMatrix(ls_pNode)
#' }

xPierMatrix <- function(list_pNode, displayBy=c("score","rank","weight","pvalue"), combineBy=c('intersect','union'), aggregateBy=c("none","fishers","logistic","Ztransform","orderStatistic"), verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    displayBy <- match.arg(displayBy)
    combineBy <- match.arg(combineBy)
    aggregateBy <- match.arg(aggregateBy) 
    
   	if(any(class(list_pNode) %in% c("pNode"))){
		list_pNode <- list(list_pNode)
	}else if(class(list_pNode)=="list"){
		## Remove null elements in a list
		list_pNode <- base::Filter(base::Negate(is.null), list_pNode)
		if(length(list_pNode)==0){
			return(NULL)
		}
	}else{
		stop("The function must apply to 'list' of 'pNode' objects or a 'pNode' object.\n")
	}
	
	## get nodes involved
	ls_nodes <- lapply(list_pNode, function(x){
		V(x$g)$name
	})
	if(combineBy=='intersect'){
		nodes <- base::Reduce(intersect, ls_nodes)
	}else if(combineBy=='union'){
		nodes <- base::Reduce(union, ls_nodes)
	}
	nodes <- sort(nodes)
	
	## Combine into a data frame called 'df_predictor'
	list_names <- names(list_pNode)
	if(is.null(list_names)){
		list_names <- paste('Predictor', 1:length(list_pNode), sep=' ')
		names(list_pNode) <- list_names
	}
	ls_priority <- lapply(list_pNode, function(pNode){
		p <- pNode$priority
		ind <- match(nodes, rownames(p))
		#ind <- ind[!is.na(ind)]
		if(displayBy=='score' | displayBy=='pvalue'){
			res <- p[ind, c("priority")]
		}else if(displayBy=='rank'){
			res <- p[ind, c("rank")]
		}else if(displayBy=='weight'){
			res <- p[ind, c("weight")]
		}
	})
	df_predictor <- do.call(cbind, ls_priority)
	rownames(df_predictor) <- nodes
	
	## replace NA with worst value
	if(displayBy=='score' | displayBy=='weight' | displayBy=='pvalue'){
		df_predictor[is.na(df_predictor)] <- 0
	}else if(displayBy=='rank'){
		df_predictor[is.na(df_predictor)] <- length(nodes)
	}
	
	## only when displayBy=='pvalue'
	## Convert into p-values by computing an empirical cumulative distribution function
	if(displayBy=='pvalue'){
		ls_pval <- lapply(1:ncol(df_predictor), function(j){
			x <- df_predictor[,j]
			my.CDF <- stats::ecdf(x)
			pval <- 1 - my.CDF(x)
		})
		df_pval <- do.call(cbind, ls_pval)
		rownames(df_pval) <- rownames(df_predictor)
		colnames(df_pval) <- colnames(df_predictor)
		df_predictor <- df_pval
		
		## aggregate p values
		if(aggregateBy != "none"){
			df_ap <- dnet::dPvalAggregate(pmatrix=df_predictor, method=aggregateBy)
			df_ap <- sort(df_ap, decreasing=FALSE)
			
			## get rank
			df_rank <- rank(df_ap, ties.method="min")
			######
			df_ap[df_ap==0] <- min(df_ap[df_ap!=0])
			######
			## adjp
			df_adjp <- stats::p.adjust(df_ap, method="BH")
			######
			## priority: first log10-transformed ap and then being rescaled into the [0,5] range
			priority <- -log10(df_ap)
			priority <- 5 * (priority - min(priority))/(max(priority) - min(priority))
			
			## df_priority
			df_priority <- data.frame(name=names(df_ap), rank=df_rank, pvalue=df_ap, fdr=df_adjp, priority=priority, stringsAsFactors=FALSE)
			### add description
			df_priority$description <- xSymbol2GeneID(df_priority$name, details=TRUE, RData.location=RData.location)$description
			###
			
			## df_predictor
			ind <- match(names(df_ap), rownames(df_predictor))
			df_predictor <- df_predictor[ind,]
			
			############
			## seed info
			predictor_names <- names(list_pNode)
			predictor_names <- gsub('^Annotation_', '', predictor_names)
			predictor_names <- gsub('_.*', '', predictor_names)
			ls_df <- lapply(1:length(list_pNode), function(i){
				pNode <- list_pNode[[i]]
				genes <- rownames(pNode$priority)[pNode$priority$seed==1]
				df <- data.frame(Gene=genes, Predictor=rep(predictor_names[i], length(genes)), stringsAsFactors=FALSE)
			})
			df <- do.call(rbind, ls_df)
			mat <- as.matrix(xSparseMatrix(df, rows=df_priority$name, columns=NULL, verbose=FALSE))
			ind_row <- match(df_priority$name, rownames(mat))
			ind_col <- match(unique(predictor_names), colnames(mat))
			mat <- mat[ind_row, ind_col]
			overall <- apply(mat!=0, 1, sum)
			############
			
			##############
			## get edges involved
			ls_edges <- lapply(list_pNode, function(x){
				relations <- igraph::get.data.frame(x$g, what="edges")
			})
			edges <- unique(do.call(rbind, ls_edges))
			## get metag
			metag <- igraph::graph.data.frame(d=edges, directed=FALSE, vertices=nodes)
			##############
			
			dTarget <- list(priority  = cbind(df_priority,Overall=overall, mat),
							predictor = df_predictor,
							metag	  = metag, 
							Call      = match.call()
						 )
			class(dTarget) <- "dTarget"
			
			df_predictor <- dTarget
		}
		
	}
	
	if(verbose){
		
		if(displayBy=="pvalue" & aggregateBy!="none"){
			message(sprintf("A total of %d genes are prioritised, combined by '%s' and aggregated by '%s' from %d predictors", nrow(df_predictor$priority), combineBy, aggregateBy, length(list_pNode)), appendLF=TRUE)
		}else{
			message(sprintf("A matrix of %d genes x %d predictors are generated, displayed by '%s' and combined by '%s'", nrow(df_predictor), ncol(df_predictor), displayBy, combineBy), appendLF=TRUE)
		}
		
	}
	
	
    invisible(df_predictor)
}
