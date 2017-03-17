#' Function to visualise evidence for prioritised genes in a gene network
#'
#' \code{xVisEvidence} is supposed to visualise evidence for prioritised genes in a gene network. It returns an object of class "igraph". 
#'
#' @param dTarget an object of class "dTarget"
#' @param g an object of class "igraph". If NA, the 'metag' will be used, which is part of the object of class "dTarget"
#' @param nodes which node genes are in query. If NULL, the top gene will be queried
#' @param neighbor.order an integer giving the order of the neighborhood. By default, it is 1-order neighborhood
#' @param neighbor.seed logical to indicate whether neighbors are seeds only. By default, it sets to true
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{http://html-color-codes.info/color-names}
#' @param legend.position the legend position. If NA, the legend will be hiden
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param ... additional graphic parameters. See \url{http://igraph.org/r/doc/plot.common.html} for the complete list.
#' @return
#' a subgraph, an object of class "igraph".
#' @export
#' @seealso \code{\link{xPierMatrix}}
#' @include xVisEvidence.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' subg <- xVisEvidence(dTarget, nodes="B2M", neighbor.order=1, neighbor.seed=TRUE, vertex.label.color="black", vertex.label.cex=1, vertex.label.dist=1, legend.position="bottomleft", signature=FALSE)
#' }

xVisEvidence <- function(dTarget, g=NA, nodes=NULL, neighbor.order=1, neighbor.seed=TRUE, colormap="ggplot2", legend.position="topleft", verbose=TRUE, ...)
{

    if(class(dTarget) == "dTarget"){
        df_evidence <- dTarget$priority[, 6:ncol(dTarget$priority)]
        
		if(!is.na(g)){
			g <- match.arg(g)
		}else{
			g <- dTarget$metag
		}
		
    }else{
    	stop("The function must apply to a 'dTarget' object.\n")
    }
	
	if(is.null(nodes)){
		nodes <- rownames(df_evidence)[1]
	}else{
		ind <- match(nodes, rownames(df_evidence))
		ind <- which(!is.na(ind))
		if(length(ind)>=1){
			nodes <- nodes[ind]
		}else{
			nodes <- rownames(df_evidence)[1]
		}
	}
    
    neighs.out <- igraph::neighborhood(g, order=neighbor.order, nodes=nodes, mode="all")
	neighbors <- names(unlist(neighs.out))
	if(neighbor.seed){
		# restrict to seeds
		ind <- neighbors %in% rownames(df_evidence)[df_evidence$Overall>0]
		neighbors <- neighbors[ind]
	}
	vids <- union(neighbors, nodes)
	#subg <- igraph::induced_subgraph(g, vids=vids, delete.vertices=TRUE)
    subg <- dnet::dNetInduce(g, nodes_query=vids, knn=0, largest.comp=TRUE)
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("The %d-order graph induced by %d input genes: %d nodes and %d edges", neighbor.order, length(nodes), vcount(subg), ecount(subg)), appendLF=TRUE)
	}
	
	ind <- match(V(subg)$name, rownames(df_evidence))
	df_val <- df_evidence[ind,]
	ls_val <- lapply(1:nrow(df_val), function(i){
		if(sum(df_val[i,-1])==0){
			NULL
		}else{
			as.numeric(df_val[i,-1]!=0)
		}
	})
	names(ls_val) <- rownames(df_val)
	
	## nodes NULL are drawn as circles
	vertex.shape <- rep("pie", length(ls_val))
	vertex.shape[sapply(ls_val, is.null)] <- "circle"
	## pie color
	pie.color <- xColormap("ggplot2")(ncol(df_val)-1)
	## legend text
	legend.text <- colnames(df_val)[-1]
	legend.text[legend.text=="OMIM"] <- "dGene"
	legend.text[legend.text=="Phenotype"] <- "pGene"
	legend.text[legend.text=="Function"] <- "fGene"
	legend.text[legend.text=="nearbyGenes"] <- "nGene"
	legend.text[legend.text=="eQTL"] <- "eGene"
	legend.text[legend.text=="HiC"] <- "hGene"
	## vertex size
	vertex.size <- igraph::degree(subg)
	if(min(vertex.size) == max(vertex.size)){
		vertex.size <- 12
	}else{
		vertex.size <- 12 * (vertex.size - min(vertex.size))/(max(vertex.size) - min(vertex.size)) + 8
	}
	## draw graph
	xVisNet(subg, vertex.shape=vertex.shape, vertex.pie=ls_val, vertex.pie.color=list(pie.color), vertex.pie.border="grey", vertex.color="grey", vertex.size=vertex.size, ...)
	if(!is.na(legend.position)){
		legend(legend.position, legend=legend.text, col=pie.color, pch=17, bty="n", pt.cex=1.2, cex=1, text.col="black", text.font=4, horiz=FALSE)
	}
	
    return(subg)
}
