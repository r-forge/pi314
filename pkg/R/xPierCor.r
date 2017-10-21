#' Function to calculate correlation between prioritised genes and user-defined external data
#'
#' \code{xPierCor} is supposed to calculate correlation between prioritised genes and user-defined external data. 
#
#' @param pNode an object of class "pNode" (or "sTarget" or "dTarget"). Alternatively, it can be a data frame with two columns ('name' and 'priority')
#' @param data a named input vector containing numeric data for genes (gene symbols)
#' @param method the method used to calcualte correlation. It can be 'pearson' for Pearson's correlation or 'spearman' for Spearman rank correlation
#' @param pvalue.type the type of the p-value calcualted. It can be 'nominal' for nominal p-value or 'empirical' for empirical p-value
#' @param seed an integer specifying the seed
#' @param nperm the number of random permutations
#' @param plot logical to indicate whether scatter plot is drawn
#' @return an object of class "ggplot" with a component 'res' (a data frame with two columns 'cor' and 'pval') if the scatter plot is drawn. Otherwise it is a data frame with two columns 'cor' and 'pval'.
#' @note none
#' @export
#' @seealso \code{\link{xPierCor}}
#' @include xPierCor.r
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
#' # c) do correlation
#' data <- pNode$priority$priority[1:100]
#' name(data) <- pNode$priority$name[1:100]
#' res <- xPierCor(pNode, data, method="pearson", pvalue.type="empirical", nperm=2000, plot=TRUE)
#' }

xPierCor <- function(pNode, data, method=c("pearson","spearman"), pvalue.type=c("nominal","empirical"), seed=825, nperm=2000, plot=FALSE)
{
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    method <- match.arg(method)
    pvalue.type <- match.arg(pvalue.type)
    
    if(class(pNode) == "pNode"){
        df_priority <- pNode$priority[, c("name","priority")]
    }else if(class(pNode) == "sTarget" | class(pNode) == "dTarget"){
    	df_priority <- pNode$priority[, c("name","priority")]
    }else if(class(pNode) == "data.frame"){
    	df_priority <- pNode[,c(1:2)]
    	colnames(df_priority) <- c("name","priority")
    }else{
    	stop("The function must apply to a 'pNode' or 'sTarget' or 'dTarget' object.\n")
    }
    
    ############
    if(length(data)==0){
    	return(NULL)
    }
    ############
    if (is.vector(data)){
    	# assume a vector
		if(is.null(names(data))){
			stop("The input data must have names with attached gene symbols.\n")
		}
    }else{
        stop("The input data must be a vector.\n")
    }
    
	df_priority_data <- df_priority
	ind <- match(df_priority_data$name, names(data))
	df_priority_data$data <- data[ind]
	df <- subset(df_priority_data, !is.na(data))
    
    ##############
    res <- stats::cor.test(x=df$priority, y=as.numeric(df$data), method=method, exact=FALSE)
	cor_obs <- signif(res$estimate, 3)
    ##############
    
    if(pvalue.type == 'nominal'){
    	pval_obs <- res$p.value
    }else if(pvalue.type == 'empirical'){
		B <- nperm
		set.seed(seed)
		vec_p <- sapply(1:B, function(i){
			df$priority <- sample(df_priority_data$priority, nrow(df))
			cor_exp <- stats::cor(x=df$priority, y=df$data, method=method)
		})
		pval_obs <- sum(abs(vec_p) > abs(cor_obs))/B
    }
    
    if(pval_obs < 0.05){
    	pval_obs <- format(signif(res$p.value,2), scientific=TRUE)
    }else{
    	pval_obs <- signif(res$p.value,3)
    }
    
    if(plot){
		priority <- name <- NULL
		m <- ggplot(df, aes(x=priority, y=data))
		m <- m + geom_point() 
		m <- m + geom_smooth(method=c("lm","loess")[1], se=TRUE, span=4)
		m <- m + theme_bw() + theme(legend.position="top", axis.title.y=element_text(size=12,color="black"), axis.text.y=element_text(size=8,color="black"), axis.title.x=element_text(size=12,color="black"), axis.text.x=element_text(size=8,color="black"), panel.background=element_rect(fill="transparent"))
		m <- m + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
		gp_curve <- m + labs(x=paste0("Pi's 5-star rating"), title=paste0("Correlation (",method,")"), subtitle=paste0("correlation: ",cor_obs,', ',pvalue.type,' p-value: ',pval_obs)) + theme(plot.title=element_text(hjust=0.5, size=12), plot.subtitle=element_text(hjust=0.5, size=10))
		
		if(0){
		gp_curve <- gp_curve + scale_x_continuous(limits=c(0,ceiling(max(df$priority)*10)/10)) 
		gp_curve <- gp_curve + scale_y_reverse(limits=c(ceiling(max(df$data)*10)/10, floor(min(df$data)*10)/10))
		gp_curve <- gp_curve + ggrepel::geom_text_repel(aes(label=name), size=2, fontface='bold', box.padding=unit(0.2,"lines"), point.padding=unit(0.2,"lines"), segment.color='grey50', segment.alpha=0.5, arrow=arrow(length=unit(0.01,'npc')), force=0.1)
		}
		
		gp_curve$res <- data.frame(cor=cor_obs, pval=pval_obs, stringsAsFactors=FALSE)
	}else{
		gp_curve <- data.frame(cor=cor_obs, pval=pval_obs, stringsAsFactors=FALSE)
	}
    
    invisible(gp_curve)
}
