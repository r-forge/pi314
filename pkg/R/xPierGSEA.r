#' Function to prioritise pathways based on GSEA analysis of prioritised genes
#'
#' \code{xPierGSEA} is supposed to prioritise pathways given prioritised genes and the ontology in query. It is done via gene set enrichment analysis (GSEA). It returns an object of class "eGSEA". 
#
#' @param pNode an object of class "pNode" (or "pTarget" or "dTarget")
#' @param priority.top the number of the top targets used for GSEA. By default, it is NULL meaning all targets are used
#' @param ontology the ontology supported currently. It can be "GOBP" for Gene Ontology Biological Process, "GOMF" for Gene Ontology Molecular Function, "GOCC" for Gene Ontology Cellular Component, "PS" for phylostratific age information, "PS2" for the collapsed PS version (inferred ancestors being collapsed into one with the known taxonomy information), "SF" for domain superfamily assignments, "DO" for Disease Ontology, "HPPA" for Human Phenotype Phenotypic Abnormality, "HPMI" for Human Phenotype Mode of Inheritance, "HPCM" for Human Phenotype Clinical Modifier, "HPMA" for Human Phenotype Mortality Aging, "MP" for Mammalian Phenotype, and Drug-Gene Interaction database (DGIdb) for drugable categories, and the molecular signatures database (Msigdb, including "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7")
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 500
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths)
#' @param weight an integer specifying score weight. It can be "0" for unweighted (an equivalent to Kolmogorov-Smirnov, only considering the rank), "1" for weighted by input gene score (by default), and "2" for over-weighted, and so on
#' @param nperm the number of random permutations. For each permutation, gene-score associations will be permutated so that permutation of gene-term associations is realised
#' @param fast logical to indicate whether to fast calculate GSEA resulting. By default, it sets to true, but not necessarily does so. It will depend on whether the package "fgsea" has been installed. It can be installed via: \code{source("http://bioconductor.org/biocLite.R"); biocLite(c("fgsea"))}. If not yet installed, this option will be disabled
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to false for no display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "eGSEA", a list with following components:
#' \itemize{
#'  \item{\code{df_summary}: a data frame of nTerm x 9 containing gene set enrichment analysis result, where nTerm is the number of terms/genesets, and the 9 columns are "setID" (i.e. "Term ID"), "name" (i.e. "Term Name"), "nAnno" (i.e. number in members annotated by a term), "nLead" (i.e. number in members as leading genes), "es" (i.e. enrichment score), "nes" (i.e. normalised enrichment score; enrichment score but after being normalised by gene set size), "pvalue" (i.e. nominal p value), "adjp" (i.e. adjusted p value; p value but after being adjusted for multiple comparisons), "distance" (i.e. term distance or metadata)}
#'  \item{\code{leading}: a list of gene sets, each storing leading gene info (i.e. the named vector with names for gene symbols and elements for priority rank). Always, gene sets are identified by "setID"}
#'  \item{\code{leading}: a list of gene sets, each storing full info on gene set enrichment analysis result (i.e. a data frame of nGene x 5, where nGene is the number of genes, and the 5 columns are "GeneID", "Rank" for priority rank, "Score" for priority score, "RES" for running enrichment score, and "Hits" for gene set hits info with 1 for gene hit, 2 for leading gene hit, 3 for the point defining leading genes, 0 for no hit). Always, gene sets are identified by "setID"}
#' }
#' @note none
#' @export
#' @seealso \code{\link{xGSEAbarplot}}, \code{\link{xGSEAdotplot}}
#' @include xPierGSEA.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#'
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' # a) provide the seed nodes/genes with the weight info
#' ## load ImmunoBase
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase', RData.location=RData.location)
#' ## get genes within 500kb away from AS GWAS lead SNPs
#' seeds.genes <- ImmunoBase$AS$genes_variants
#' ## seeds weighted according to distance away from lead SNPs
#' data <- 1- seeds.genes/500000
#'
#' # b) perform priority analysis
#' pNode <- xPierGenes(data=data, network="PCommonsDN_medium",restart=0.7, RData.location=RData.location)
#' 
#' # c) do pathway-level priority using GSEA
#' eGSEA <- xPierGSEA(pNode=pNode, ontology="DGIdb", nperm=2000, RData.location=RData.location)
#' bp <- xGSEAbarplot(eGSEA, top_num="auto", displayBy="nes")
#' gp <- xGSEAdotplot(eGSEA, top=1)
#' }

xPierGSEA <- function(pNode, priority.top=NULL, ontology=c("GOBP","GOMF","GOCC","PS","PS2","SF","Pfam","DO","HPPA","HPMI","HPCM","HPMA","MP", "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP", "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR", "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC", "MsigdbC6", "MsigdbC7", "DGIdb", "GTExV4", "GTExV6", "CreedsDisease", "CreedsDiseaseUP", "CreedsDiseaseDN", "CreedsDrug", "CreedsDrugUP", "CreedsDrugDN", "CreedsGene", "CreedsGeneUP", "CreedsGeneDN"), size.range=c(10,500), path.mode=c("all_paths","shortest_paths","all_shortest_paths"), weight=1, nperm=2000, fast=TRUE, verbose=TRUE, RData.location="http://galahad.well.ox.ac.uk/bigdata")
{
    startT <- Sys.time()
    message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    message("", appendLF=TRUE)
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    ontology <- match.arg(ontology)
    path.mode <- match.arg(path.mode)
    
    weight <- as.integer(weight)
    
    if(class(pNode) == "pNode"){
        df_priority <- pNode$priority[, c("seed","weight","priority","rank")]
    }else if(class(pNode) == "pTarget" | class(pNode) == "dTarget"){
    	df_priority <- pNode$priority[, c("pvalue","fdr","priority","rank")]
    }else{
    	stop("The function must apply to a 'pNode' or 'pTarget' or 'dTarget' object.\n")
    }
    
    ###############
	## priority top
	if(!is.null(priority.top)){
		priority.top <- as.integer(priority.top)
		if(priority.top > nrow(df_priority)){
			priority.top <- nrow(df_priority)
		}else if(priority.top <= 1){
			priority.top <- nrow(df_priority)
		}
	}else{
		priority.top <- nrow(df_priority)
	}
    df_priority <- df_priority[1:priority.top,]
    ###############
        
    ## convert gene symbols to entrez geneid
	name_GeneID <- xSymbol2GeneID(rownames(df_priority), check.symbol.identity=FALSE, verbose=verbose, RData.location=RData.location)
	
	ind <- which(!is.na(name_GeneID))
    data <- data.frame(priority=df_priority$priority[ind], row.names=name_GeneID[ind], stringsAsFactors=FALSE)
    
    #############################################################################################

    if(!is.na(ontology)){
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("Load the ontology %s and its gene annotations (%s) ...", ontology, as.character(now)), appendLF=TRUE)
		}

		#########
		## load GS information
		## flag the simplified version of PS
		flag_PS2 <- FALSE
		if(ontology=="PS2"){
			flag_PS2 <- TRUE
			ontology <- "PS"
		}
		GS <- xRDataLoader(RData.customised=paste('org.Hs.eg', ontology, sep=''), RData.location=RData.location, verbose=verbose)
		
		################
		if(flag_PS2){
			tmp <- as.character(unique(GS$set_info$name))
			inds <- sapply(tmp,function(x) which(GS$set_info$name==x))
		
			## new set_info
			set_info <- data.frame()
			for(i in 1:length(inds)){
				set_info<- rbind(set_info,as.matrix(GS$set_info[max(inds[[i]]),]))
			}
			## new gs
			gs <- list()
			for(i in 1:length(inds)){
				gs[[i]] <- unlist(GS$gs[inds[[i]]], use.names=FALSE)
			}
			names(gs) <- rownames(set_info)
		
			## new GS
			GS$set_info <- set_info
			GS$gs <- gs
		}
		################
		
		#########
		## get annotation information
		anno <- GS$gs

		#########
		## get ontology information
		## check the eligibility for the ontology
		all.ontologies <- c("GOBP","GOMF","GOCC","DO","HPPA","HPMI","HPCM","HPMA","MP")
		flag_ontology <- ontology %in% all.ontologies
    	
    	if(flag_ontology){
			g <- xRDataLoader(RData.customised=paste('ig.', ontology, sep=''), RData.location=RData.location, verbose=verbose)
			true.path.rule <- TRUE
			
		}else{
			nodes <- data.frame(name=as.character(GS$set_info$setID), term_id=as.character(GS$set_info$setID), term_name=as.character(GS$set_info$name), term_distance=as.character(GS$set_info$distance), stringsAsFactors=FALSE)
			nodes <- rbind(nodes, c('root','root','root','root'))
			relations <- data.frame(from='root', to=nodes$name)
			g <- igraph::graph.data.frame(d=relations, directed=TRUE, vertices=nodes)
			true.path.rule <- FALSE
		}
	
		# obtain the induced subgraph according to the input annotation data based on shortest paths (i.e. the most concise subgraph induced)
		subg <- xDAGanno(g=g, annotation=anno, path.mode=path.mode, true.path.rule=true.path.rule, verbose=verbose)
		anno <- V(subg)$anno
		names(anno) <- V(subg)$name
		g <- subg
	
	}else{
		stop("There is no input for the ontology.\n")
	}
    #############################################################################################
    
    flag_fgsea <- FALSE
    pkgs <- c("fgsea")
    if(all(pkgs %in% rownames(utils::installed.packages()))){
        tmp <- sapply(pkgs, function(pkg) {
            requireNamespace(pkg, quietly=TRUE)
        })
        if(all(tmp)){
        	flag_fgsea <- TRUE
        }
    }
    
    if(flag_fgsea & fast){
    
		if(verbose){
			now <- Sys.time()
			message(sprintf("\n#######################################################", appendLF=TRUE))
			message(sprintf("'fgsea' from the fgsea package is being called (%s):", as.character(now)), appendLF=TRUE)
			message(sprintf("#######################################################", appendLF=TRUE))
		}
		
		stats <- data$priority
		names(stats) <- rownames(data)
		fgseaRes <- fgsea::fgsea(pathway=anno, stats=stats, minSize=size.range[1], maxSize=size.range[2], gseaParam=weight, nperm=nperm)
		tab <- data.frame(setID         = fgseaRes$pathway,
						   ES           = fgseaRes$ES,
						   nES          = fgseaRes$NES,
						   pvalue       = fgseaRes$pval,
						   adjp         = fgseaRes$padj,
						   setSize      = fgseaRes$size,
						   stringsAsFactors=FALSE
						  )
		pvalue <- nES <- NULL
		res <- tab[with(tab,order(pvalue,-nES)),]
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("#######################################################", appendLF=TRUE))
			message(sprintf("'fgsea' has been finished (%s)!", as.character(now)), appendLF=TRUE)
			message(sprintf("#######################################################\n", appendLF=TRUE))
		}
	
	}else{
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("\n#######################################################", appendLF=TRUE))
			message(sprintf("'dGSEA' from the dnet package is being called (%s):", as.character(now)), appendLF=TRUE)
			message(sprintf("#######################################################", appendLF=TRUE))
		}
	
		eTerm <- dGSEA(data=data, identity="entrez", check.symbol.identity=FALSE, ontology="Customised", customised.genesets=anno, sizeRange=size.range, which_distance=NULL, weight=weight, nperm=nperm, fast=TRUE, sigTail="one-tail", p.adjust.method="BH", verbose=verbose, RData.location=RData.location)
	
		res <- dGSEAview(eTerm, which_sample=1, top_num=NULL, sortBy="pvalue", decreasing=TRUE, details=FALSE)
		res <- res[,c("setID","ES","nES","pvalue","adjp","setSize")]
		rownames(res) <- NULL
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("#######################################################", appendLF=TRUE))
			message(sprintf("'dGSEA' has been finished (%s)!", as.character(now)), appendLF=TRUE)
			message(sprintf("#######################################################\n", appendLF=TRUE))
		}
	}

    #############################################################################################
	### get Leading genes
    if(TRUE){

		geneid <- rownames(data)
		nGene <- nrow(data)
		## score rank
		rank.score <- data$priority
		ind <- order(rank.score, decreasing=TRUE)
		rank.score.sorted <- rank.score[ind]
		geneid.sorted <- geneid[ind]
    
		ls_df_leading <- lapply(res$setID, function(x){
			#df_leading <- visGSEA(eTerm, which_sample=1, which_term=x, plot=FALSE)
			
			## initialisation
			nHit <- length(anno[[x]])
			nMiss <- nGene - nHit
			## observed
			observed.point <- rep(-1/nMiss, nGene)
			flag <- match(anno[[x]], geneid.sorted)
			###### remove NA
			flag <- flag[!is.na(flag)]
			######
			if(weight==0) {
				observed.point[flag] <- 1/nHit
			}else if(weight==1){
				hit_tmp <- abs(rank.score.sorted[flag])
				observed.point[flag] <- hit_tmp/sum(hit_tmp)
			}else{
				hit_tmp <- abs(rank.score.sorted[flag] ** weight)
				observed.point[flag] <- hit_tmp/sum(hit_tmp)      
			}
			RES <- cumsum(observed.point)
			max.RES <- max(RES)
			min.RES <- min(RES)
			es.observed <- signif(ifelse(max.RES>abs(min.RES), max.RES, min.RES), digits=5)
			es.position <- ifelse(max.RES>abs(min.RES), which.max(RES), which.min(RES))
			## for leading genes
			if(RES[es.position]<0){
				ind <- which(flag >= es.position)
			}else{
				ind <- which(flag <= es.position)
			}
			hits <- rep(0, length(RES))
			hits[flag] <- 1
			hits[flag[ind]] <- 2
			hits[es.position] <- 3
			df_leading <- data.frame(GeneID=geneid.sorted, Rank=1:length(RES), Score=rank.score.sorted, RES=RES, Hits=hits, stringsAsFactors=FALSE)
			
		})
		names(ls_df_leading) <- res$setID
		
		ls_leadingGenes <- lapply(ls_df_leading, function(x){
			x$GeneID[x$Hits==2]
		})
		
    }
	
	# replace EntrezGenes with gene symbols	
	if(1){
		## load Enterz Gene information
		EG <- xRDataLoader(RData.customised=paste('org.Hs.eg', sep=''), RData.location=RData.location, verbose=verbose)
		allGeneID <- EG$gene_info$GeneID
		allSymbol <- as.vector(EG$gene_info$Symbol)
		
		leadingGenes <- lapply(ls_leadingGenes,function(x){
			ind <- match(x, allGeneID)
			y <- allSymbol[ind]
			
			ind <- match(y, rownames(df_priority))
			rank <- df_priority$rank[ind]
			names(rank) <- y
			return(rank)
		})
		
		## append leading genes
		res$nLead <- sapply(leadingGenes,length)
	}
	
	## append "term_name" and "term_distance"
	ind <- match(res$setID, V(g)$name)
	summary <- data.frame(setID=res$setID, name=V(g)$term_name[ind], nAnno=res$setSize, nLead=res$nLead, es=res$ES, nes=res$nES, pvalue=res$pvalue, adjp=res$adjp, distance=V(g)$term_distance[ind], stringsAsFactors=FALSE)
	
	## scientific notation
	summary$es <- signif(summary$es, digits=3)
	summary$nes <- signif(summary$nes, digits=3)
	summary$pvalue <- signif(summary$pvalue, digits=3)
	summary$pvalue <- ifelse(summary$pvalue<0.01 & summary$pvalue!=0, as.numeric(format(summary$pvalue,scientific=TRUE)), summary$pvalue)
	summary$adjp <- signif(summary$adjp, digits=3)
	summary$adjp <- ifelse(summary$adjp<0.01 & summary$adjp!=0, as.numeric(format(summary$adjp,scientific=TRUE)), summary$adjp)
	
    eGSEA <- list(df_summary = summary,
    			  leading = leadingGenes,
    			  full = ls_df_leading,
    			  Call = match.call()
                 )
    class(eGSEA) <- "eGSEA"
    
####################################################################################
    endT <- Sys.time()
    message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    
    invisible(eGSEA)
}
