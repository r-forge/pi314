#' Function to simulate gold standard negatives (GSN) given gold standard positives (GSP) and a gene network
#'
#' \code{xGSsimulator} is supposed to simulate gold standard negatives (GSN) given gold standard positives (GSP) and an input gene network. GSN targets are those after excluding GSP targets and their 1-order (by default) neighbors in the gene network.
#'
#' @param GSP a vector containing Gold Standard Positives (GSP)
#' @param population a vector containing population space in which gold standard negatives (GSN) will be considered. By default, it is NULL, meaning genes in the network will be used instead
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathways Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). By default, "STRING_medium" and "PCommonsUN_medium" are used
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network
#' @param neighbor.order an integer giving the order of the neighborhood. By default, it is 1-order neighborhood
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' a list with following components:
#' \itemize{
#'  \item{\code{GSN}: a vector containing simulated GSN}
#'  \item{\code{GSP}: a vector containing GSP after considering the population space}
#'  \item{\code{g}: an "igraph" object}
#' }
#' @note If multiple graphs are provided, they will be unionised for use.
#' @export
#' @seealso \code{\link{xRDataLoader}}, \code{\link{xPredictROCR}}, \code{\link{xMLrandomforest}}
#' @include xGSsimulator.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' sGS <- xGSsimulator(GSP, population, network=c("STRING_medium","PCommonsUN_medium"), RData.location=RData.location)
#' }

xGSsimulator <- function(GSP, population=NULL, network=c("STRING_medium","STRING_low","STRING_high","STRING_highest","PCommonsUN_high","PCommonsUN_medium")[c(1,6)], network.customised=NULL, neighbor.order=1, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{

    if(is.null(GSP)){
        stop("The input GSP must be not NULL.\n")
    }
    
    if(!is.null(network.customised) && class(network.customised)=="igraph"){
		if(verbose){
			now <- Sys.time()
			message(sprintf("Load the customised network (%s) ...", as.character(now)), appendLF=TRUE)
		}
		ig <- network.customised
		
	}else{
	
		default.network <- c("STRING_medium","STRING_low","STRING_high","STRING_highest","PCommonsUN_high","PCommonsUN_medium")
		ind <- match(default.network, network)
		network <- default.network[!is.na(ind)]
	
		if(length(network) > 0){
			###########################	
			# built-in network
			###########################
			if(length(grep('STRING',network,perl=TRUE)) > 0){
				STRING <- xRDataLoader(RData.customised='org.Hs.string', verbose=verbose, RData.location=RData.location)
			}
			if(length(grep('PCommonsUN',network,perl=TRUE)) > 0){
				PC <- xRDataLoader(RData.customised='org.Hs.PCommons_UN', verbose=verbose, RData.location=RData.location)
			}
			
			ig_list <- lapply(network, function(x){
				if(verbose){
					now <- Sys.time()
					message(sprintf("Load the network %s (%s) ...", x, as.character(now)), appendLF=TRUE)
				}
				
				if(length(grep('STRING',x,perl=TRUE)) > 0){
					g <- STRING
			
					## restrict to those edges with given confidence
					flag <- unlist(strsplit(x,"_"))[2]
					if(flag=='highest'){
						eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=900])"))
					}else if(flag=='high'){
						eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=700])"))
					}else if(flag=='medium'){
						eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=400])"))
					}else if(flag=='low'){
						eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[combined_score>=150])"))
					}
			
					########################
					# because of the way storing the network from the STRING database
					## extract relations (by symbol)
					V(g)$name <- V(g)$symbol
					relations <- igraph::get.data.frame(g, what="edges")[, c(1,2)]
					relations <- unique(relations)
					colnames(relations) <- c("from","to")
					########################
			
					g <- igraph::graph.data.frame(d=relations, directed=FALSE)
			
				}else if(length(grep('PCommonsUN',x,perl=TRUE)) > 0){
					g <- PC
			
					flag <- unlist(strsplit(x,"_"))[2]
					if(flag=='high'){
						# restrict to those edges with physical interactions and with score>=102
						eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[in_complex_with>=102 | interacts_with>=102])"))
					}else if(flag=='medium'){
						# restrict to those edges with physical interactions and with score>=101
						eval(parse(text="g <- igraph::subgraph.edges(g, eids=E(g)[in_complex_with>=101 | interacts_with>=101])"))
					}
			
					relations <- igraph::get.data.frame(g, what="edges")[, c(1,2)]
					relations <- unique(relations)
					colnames(relations) <- c("from","to")
					g <- igraph::graph.data.frame(d=relations, directed=FALSE)
			
				}
				
				return(g)	
			})
			# unionise graphs
			ig <- igraph::graph.union(ig_list)
		
		}else{
			ig <- NULL
		}
	
	}
    ######################################################################################
    
    if(is.null(ig)){
        stop("The input network must be not NULL.\n")
    }
    
	## neighbor to be excluded
	nodes <- GSP[which(GSP %in% V(ig)$name)]
	neighs.out <- igraph::neighborhood(ig, order=neighbor.order, nodes=nodes, mode="all")
	neigh_excluded <- unique(V(ig)[unlist(neighs.out)]$name)
    
    ## consider whether the population space is provided
    if(is.null(population)){
    	population <- base::Reduce(union, list(V(ig)$name,GSP))
    }
    
    ## GSP
    res_GSP <- base::Reduce(intersect, list(population,GSP))
    
    ## GSP neighbors
    res_neigh_excluded <- base::Reduce(intersect, list(population,neigh_excluded))
    
    ## GSN
    excluded <- base::Reduce(union, list(res_neigh_excluded,res_GSP))
    GSN <- base::setdiff(population, excluded)
    
	if(verbose){
		message(sprintf("Amongst %d genes from populations, %d genes (out of %d input genes) in GSP, %d genes as %d-order GSP neighbors, %d genes in GSN", length(population), length(res_GSP), length(GSP), length(res_neigh_excluded), neighbor.order, length(GSN)), appendLF=TRUE)
	}
    
    sGS <- list(
    			GSN = GSN,
    			GSP = res_GSP,
               	g = ig
              )
       
    invisible(sGS)
}
