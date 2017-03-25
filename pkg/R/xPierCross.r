#' Function to extract priority matrix from a list of dTarget objects
#'
#' \code{xPierCross} is supposed to extract priority matrix from a list of dTarget objects. Also supported is the aggregation of priority matrix (similar to the meta-analysis) generating the priority results; we view this functionality as the cross mode of the prioritisation.
#'
#' @param list_dTarget a list of "dTarget" objects or a "dTarget" object
#' @param displayBy which priority will be extracted. It can be "priority" for priority score (by default), "rank" for priority rank, "pvalue" for priority p-value, "fdr" for priority fdr
#' @param combineBy how to resolve nodes/targets from a list of "dTarget" objects. It can be "intersect" for intersecting nodes (by default), "union" for unionising nodes
#' @param aggregateBy the aggregate method used. It can be either "none" for no aggregation, or "orderStatistic" for the method based on the order statistics of p-values, "fishers" for Fisher's method, "Ztransform" for Z-transform method, "logistic" for the logistic method. Without loss of generality, the Z-transform method does well in problems where evidence against the combined null is spread widely (equal footings) or when the total evidence is weak; Fisher's method does best in problems where the evidence is concentrated in a relatively small fraction of the individual tests or when the evidence is at least moderately strong; the logistic method provides a compromise between these two. Notably, the aggregate methods 'fishers' and 'logistic' are preferred here
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' If aggregateBy is 'none' (by default), a data frame containing priority matrix, with each column/disease for either priority score, or priorty rank or priority p-value.
#' If aggregateBy is not 'none', an object of the class "cTarget", a list with following components:
#' \itemize{
#'  \item{\code{priority}: a data frame of nGene X 6 containing gene priority (aggregated) information, where nGene is the number of genes, and the 6 columns are "name" (gene names), "rank" (ranks of the priority scores), "pvalue" (the aggregated p-value, converted from empirical cumulative distribution of the probability of being GSP), "fdr" (fdr adjusted from the aggregated p-value), "priority" (-log10(pvalue) but rescaled into the 5-star ratings), "description" (gene description)}
#'  \item{\code{disease}: a data frame containing disease matrix, with each column/disease for either priority score, or priorty rank or priority p-value}
#'  \item{\code{call}: the call that produced this result}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xPierMatrix}}
#' @include xPierCross.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' df_score <- xPierCross(ls_dTarget)
#' }

xPierCross <- function(list_dTarget, displayBy=c("priority","rank","pvalue","fdr"), combineBy=c('intersect','union'), aggregateBy=c("none","fishers","logistic","Ztransform","orderStatistic"), verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    displayBy <- match.arg(displayBy)
    combineBy <- match.arg(combineBy)
    aggregateBy <- match.arg(aggregateBy) 
    
   	if(any(class(list_dTarget) %in% c("dTarget","pTarget"))){
		list_dTarget <- list(list_dTarget)
	}else if(class(list_dTarget)=="list"){
		## Remove null elements in a list
		list_dTarget <- base::Filter(base::Negate(is.null), list_dTarget)
		if(length(list_dTarget)==0){
			return(NULL)
		}
	}else{
		stop("The function must apply to 'list' of 'dTarget' objects or a 'pTarget' object.\n")
	}
	
	## get nodes involved
	ls_nodes <- lapply(list_dTarget, function(x){
		x$priority$name
	})
	if(combineBy=='intersect'){
		nodes <- base::Reduce(intersect, ls_nodes)
	}else if(combineBy=='union'){
		nodes <- base::Reduce(union, ls_nodes)
	}
	nodes <- sort(nodes)
	
	## Combine into a data frame called 'df_disease'
	list_names <- names(list_dTarget)
	if(is.null(list_names)){
		list_names <- paste('Disease', 1:length(list_dTarget), sep=' ')
		names(list_dTarget) <- list_names
	}
	ls_priority <- lapply(list_dTarget, function(dTarget){
		p <- dTarget$priority
		ind <- match(nodes, rownames(p))
		res <- p[ind, displayBy]
	})
	df_disease <- do.call(cbind, ls_priority)
	rownames(df_disease) <- nodes
	
	## replace NA with worst value
	if(displayBy=='priority'){
		df_disease[is.na(df_disease)] <- 0
	}else if(displayBy=='pvalue' | displayBy=='fdr'){
		df_disease[is.na(df_disease)] <- 1
	}else if(displayBy=='rank'){
		df_disease[is.na(df_disease)] <- length(nodes)
	}
	
	## only when displayBy=='pvalue'
	## Convert into p-values by computing an empirical cumulative distribution function
	if(displayBy=='pvalue'){

		## aggregate p values
		if(aggregateBy != "none"){
			df_ap <- dnet::dPvalAggregate(pmatrix=df_disease, method=aggregateBy)
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
			df_priority$description <- XGR::xSymbol2GeneID(df_priority$name, details=TRUE, RData.location=RData.location)$description
			###
			
			## df_disease
			ind <- match(names(df_ap), rownames(df_disease))
			df_disease <- df_disease[ind,]
			
			cTarget <- list(priority = df_priority,
							predictor = df_disease, 
							Call     = match.call()
						 )
			class(dTarget) <- "cTarget"
			
			df_disease <- cTarget
		}
		
	}
	
	if(verbose){
		
		if(displayBy=="pvalue" & aggregateBy!="none"){
			message(sprintf("A total of %d genes are prioritised, combined by '%s' and aggregated by '%s' from %d predictors", nrow(df_disease$priority), combineBy, aggregateBy, length(list_dTarget)), appendLF=TRUE)
		}else{
			message(sprintf("A matrix of %d genes x %d predictors are generated, displayed by '%s' and combined by '%s'", nrow(df_disease), ncol(df_disease), displayBy, combineBy), appendLF=TRUE)
		}
		
	}
	
	
    invisible(df_disease)
}
