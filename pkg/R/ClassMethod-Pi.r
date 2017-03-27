######################################################################
#' @title Definition for S3 class \code{eTarget}
#' @description \code{eTarget} has two components: evidence and metag.
#' @param evidence a data frame
#' @param metag an 'igraph' object
#' @return an object of S3 class \code{eTarget}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' eTarget(evidence, metag)
#' }
eTarget <- function(evidence, metag){
	## integrity checks
	if(class(evidence)!='data.frame' | class(metag)!='igraph'){
		stop("The S3 class 'eTarget' object failed to pass integrity checks!\n")
	}
	value <- list(evidence=evidence, metag=metag)
	class(value) <- "eTarget"
	return(value)
}

#' @param x an object of class \code{eTarget}
#' @param ... other parameters
#' @rdname eTarget
#' @export
print.eTarget <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with two components:", class(x)), "\n", sep="")
	cat(sprintf("  $metag: an igraph object with %d nodes and %d edges", vcount(x$metag), ecount(x$metag)), "\n", sep="")
	cat(sprintf("  $evidence: a data frame of %d rows X %d columns", dim(x$evidence)[1],dim(x$evidence)[2]), "\n", sep="")
	cat("\n--------------------------------------------------\n")
	cat("$evidence:\n")
	print.data.frame(x$evidence[1:2,], row.names=FALSE)
	cat("......\n")
}


######################################################################
#' @title Definition for S3 class \code{pNode}
#' @description \code{pNode} has 7 components: priority, g, SNP, Gene2SNP, nGenes, eGenes and cGenes.
#' @param priority a data frame
#' @param g an 'igraph' object
#' @param SNP a data frame
#' @param Gene2SNP a data frame
#' @param nGenes a data frame
#' @param eGenes a data frame
#' @param cGenes a data frame
#' @return an object of S3 class \code{pNode}
#' @keywords S3 classes
#' @export
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' pNode(evidence, metag)
#' }
pNode <- function(priority, g, SNP, Gene2SNP, nGenes, eGenes, cGenes){
	## integrity checks
	if(class(priority)!='data.frame' | class(g)!='igraph'){
		stop("The S3 class 'pNode' object failed to pass integrity checks!\n")
	}
	value <- list(priority=priority, g=g, SNP=SNP, Gene2SNP=Gene2SNP, nGenes=nGenes, eGenes=eGenes, cGenes=cGenes)
	class(value) <- "pNode"
	return(value)
}

#' @param x an object of class \code{pNode}
#' @param ... other parameters
#' @rdname pNode
#' @export
print.pNode <- function(x, ...) {
	cat(sprintf("An object of S3 class '%s', with following components:", class(x)), "\n", sep="")
	cat(sprintf("  $g: an igraph object with %d nodes and %d edges", vcount(x$g), ecount(x$g)), "\n", sep="")
	cat(sprintf("  $priority: a data frame of %d rows X %d columns", dim(x$priority)[1],dim(x$priority)[2]), "\n", sep="")
	if(!is.null(x$SNP)){
		cat(sprintf("  $SNP: a data frame of %d rows X %d columns", dim(x$SNP)[1],dim(x$SNP)[2]), "\n", sep="")
	}
	if(!is.null(x$Gene2SNP)){
		cat(sprintf("  $Gene2SNP: a data frame of %d rows X %d columns", dim(x$Gene2SNP)[1],dim(x$Gene2SNP)[2]), "\n", sep="")
	}
	if(!is.null(x$nGenes)){
		cat(sprintf("  $nGenes: a data frame of %d rows X %d columns", dim(x$nGenes)[1],dim(x$nGenes)[2]), "\n", sep="")
	}
	if(!is.null(x$eGenes)){
		cat(sprintf("  $eGenes: a data frame of %d rows X %d columns", dim(x$eGenes)[1],dim(x$eGenes)[2]), "\n", sep="")
	}
	if(!is.null(x$cGenes)){
		cat(sprintf("  $cGenes: a data frame of %d rows X %d columns", dim(x$cGenes)[1],dim(x$cGenes)[2]), "\n", sep="")
	}
	cat("\n--------------------------------------------------\n")
	cat("$priority:\n")
	print.data.frame(x$priority[1:2,], row.names=FALSE)
	cat("......\n")
	if(!is.null(x$SNP)){
		cat("$SNP:\n")
		print.data.frame(x$SNP[1:2,], row.names=FALSE)
		cat("......\n")
	}
	if(!is.null(x$Gene2SNP)){
		cat("$Gene2SNP:\n")
		print.data.frame(x$Gene2SNP[1:2,], row.names=FALSE)
		cat("......\n")
	}
	if(!is.null(x$nGenes)){
		cat("$nGenes:\n")
		print.data.frame(x$nGenes[1:2,], row.names=FALSE)
		cat("......\n")
	}
	if(!is.null(x$eGenes)){
		cat("$eGenes:\n")
		print.data.frame(x$eGenes[1:2,], row.names=FALSE)
		cat("......\n")
	}
	if(!is.null(x$cGenes)){
		cat("$cGenes:\n")
		print.data.frame(x$cGenes[1:2,], row.names=FALSE)
		cat("......\n")
	}
}

