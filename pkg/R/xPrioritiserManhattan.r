#' Function to visualise prioritised genes using manhattan plot
#'
#' \code{xPrioritiserManhattan} is supposed to visualise prioritised genes using manhattan plot. Genes with the top priority are highlighed. It returns an object of class "ggplot".
#'
#' @param pNode an object of class "pNode"
#' @param color a character vector for point colors to alternate
#' @param cex a numeric value for point size
#' @param highlight.top the number of the top targets to be highlighted
#' @param highlight.col the highlight colors
#' @param highlight.label.size the highlight label size
#' @param highlight.label.offset the highlight label offset
#' @param highlight.label.col the highlight label color
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_genes', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Even the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location"
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xPrioritiser}}, \code{\link{xPrioritiserSNPs}}, \code{\link{xPrioritiserGenes}}, \code{\link{xPrioritiserPathways}}
#' @include xPrioritiserManhattan.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#'
#' # a) provide the SNPs with the significance info
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' data.file <- "http://galahad.well.ox.ac.uk/bigdata/AS.txt"
#' AS <- read.delim(data.file, header=TRUE, stringsAsFactors=FALSE)
#'
#' # b) perform priority analysis
#' pNode <- xPrioritiserSNPs(data=AS, network="PCommonsUN_medium",restart=0.7)
#'
#' # c) manhattan plot
#' mp <- xPrioritiserManhattan(pNode, highlight.top=10)
#' #pdf(file="Gene_manhattan.pdf", height=6, width=12, compress=TRUE)
#' print(mp)
#' #dev.off()
#' }

xPrioritiserManhattan <- function(pNode, color=c("darkred","darkgreen"), cex=0.5, highlight.top=20, highlight.col="deepskyblue", highlight.label.size=2, highlight.label.offset=0.02, highlight.label.col="darkblue", GR.Gene="UCSC_genes", verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/Portal")
{
    if (class(pNode) != "pNode" ){
        stop("The function must apply to a 'pNode' object.\n")
    }
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("Load positional information for Genes (%s) ...", as.character(now)), appendLF=T)
	}
    gr_Gene <- xRDataLoader(RData.customised=GR.Gene, verbose=verbose, RData.location=RData.location)
    if(is.null(gr_Gene)){
    	GR.Gene <- "UCSC_genes"
		if(verbose){
			message(sprintf("Instead, %s will be used", GR.Gene), appendLF=T)
		}
    	gr_Gene <- xRDataLoader(RData.customised=GR.Gene, verbose=verbose, RData.location=RData.location)
    }
    
    ## ONLY restricted to genes with genomic locations
	ind <- match(rownames(pNode$priority), names(gr_Gene))
	p_gr <- gr_Gene[ind[!is.na(ind)],]
	p_matrix <- pNode$priority[!is.na(ind),]
	
	## append genomic locations to GR object
	gr <- p_gr
	GenomicRanges::mcols(gr) <- cbind(GenomicRanges::mcols(gr), p_matrix[,2:4])
	## for sorting
	chrlabs <- paste('chr', as.character(c(1:22,'X','Y')), sep='')
	#eval(parse(text=paste("seqlevels(gr) <- chrlabs",sep="")))
	GenomeInfoDb::seqlevels(gr) <- chrlabs
	
	## seed only
	if(0){
		ind <- which(GenomicRanges::mcols(gr)$seed==1)
		gr <- gr[ind,]
	}
	
	## highlight points
	highlight.top <- as.integer(highlight.top)
    if ( highlight.top > length(gr) ){
        highlight.top <- length(gr)
    }
	df <- data.frame(index=1:length(gr), val=GenomicRanges::mcols(gr)$priority)
	ind_o <- df[with(df,order(-val))[1:highlight.top],1]
	gro <- gr[ind_o,]
	names(gro) <- GenomicRanges::mcols(gro)$Symbol
	
	## draw plot
    suppressWarnings(
    mp <- ggbio::plotGrandLinear(gr, eval(parse(text=paste("aes(y=priority)",sep=""))), color=color, spaceline=T, cex=cex, ylab='Priority', highlight.gr=gro, highlight.col=highlight.col, highlight.label=F, highlight.label.size=highlight.label.size, highlight.label.offset=highlight.label.offset, highlight.label.col=highlight.label.col) + theme(axis.text.x=element_text(angle=45, hjust=1,color="black",size=12), panel.background=element_rect(fill=rgb(0.95,0.95,0.95,1)))
    )
	
    invisible(mp)
}
