#' Function to visualise machine learning results
#'
#' \code{xMLplot} is supposed to visualise machine learning results. It returns an object of class "ggplot".
#'
#' @param pTarget an object of class "pTarget"
#' @param displayBy which statistics will be used for displaying. It can be either statistics across folds ("importance2fold" for predictor importance, "roc2fold" for AUC in ROC, "fmax2fold" for F-max in Precision-Recall curve) or overall statistics ("importance_accurancy" for predictor importance measured by accuracy decrease, "importance_gini" for predictor importance measured by Gini decrease, "ROC" for AUC in ROC, "Fmax" for F-max in Precision-Recall curve)
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xMLrandomforest}}
#' @include xMLplot.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
#' gp <- xMLplot(pTarget, displayBy="importance_accurancy")
#' gp
#' }

xMLplot <- function(pTarget, displayBy=c("importance2fold","roc2fold","fmax2fold","importance_accurancy","importance_gini","ROC","Fmax"), signature=TRUE) 
{
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    displayBy <- match.arg(displayBy)

    if(class(pTarget) != "pTarget"){
    	stop("The function must apply to a 'pTarget' object.\n")
    }
    
    if(displayBy == "importance2fold"){
    	df <- pTarget$importance2fold
    	xlab <- "Decrease in accuracy across folds\n(a measure of predictor importance)"
    }else if(displayBy=='roc2fold'){
    	df <- pTarget$roc2fold
    	xlab <- "AUC across folds\n(a measure of ROC)"
    }else if(displayBy=='fmax2fold'){
    	df <- pTarget$fmax2fold
    	xlab <- "F-max across folds\n(a measure of Precision-Recall curve)"
    }else if(displayBy=='importance_accurancy'){
    	df <- data.frame(Val=pTarget$importance[,1], stringsAsFactors=FALSE)
    	rownames(df) <- rownames(pTarget$importance)
    	xlab <- "Decrease in accuracy\n(a measure of predictor importance)"
    }else if(displayBy=='importance_gini'){
    	df <- data.frame(Val=pTarget$importance[,2], stringsAsFactors=FALSE)
    	rownames(df) <- rownames(pTarget$importance)
    	xlab <- "Decrease in Gini\n(a measure of predictor importance)"
    }else if(displayBy=='ROC'){
    	df <- data.frame(Val=pTarget$performance[,1], stringsAsFactors=FALSE)
    	rownames(df) <- rownames(pTarget$performance)
    	xlab <- "AUC\n(a measure of ROC)"
    }else if(displayBy=='Fmax'){
    	df <- data.frame(Val=pTarget$performance[,2], stringsAsFactors=FALSE)
    	rownames(df) <- rownames(pTarget$performance)
    	xlab <- "F-max\n(a measure of Precision-Recall curve)"
    }
    
    median <- name <- max <- min <- Val <- ''
    
    ## extract info on 'Predictor' and 'Method'
	tmp <- rownames(df)
	tmp <- gsub('^Integrated_', 'Integrated\n(', tmp)
	tmp <- gsub('^Annotation_', 'Annotation\n(', tmp)
	tmp <- gsub('^nearbyGenes_', 'nearbyGenes\n(', tmp)
	tmp <- gsub('^eQTL_', 'eQTL\n(', tmp)
	tmp <- gsub('^HiC_', 'HiC\n(', tmp)
	tmp <- gsub('^Expression_', 'Expression\n(', tmp)
	tmp <- paste(tmp,')',sep='')
	df <- data.frame(name=tmp, df, stringsAsFactors=FALSE)
	Predictor <- gsub('\n.*', '', as.character(df$name), perl=TRUE)
	Method <- gsub('.*\n\\(|\\)', '', as.character(df$name), perl=TRUE)
	df$Predictor <- Predictor
	df$Method <- Method
	
	if(displayBy == "importance2fold" | displayBy == "roc2fold" | displayBy == "fmax2fold"){
		## order by 'Predictor', 'median'
		df <- df[with(df,order(Predictor,-median)),]
		df$Predictor <- factor(df$Predictor, levels=unique(df$Predictor))
		df$Method <- factor(df$Method, levels=rev(unique(df$Method)))
		bp <- ggplot(df, aes(median, Method, colour=Predictor))
		bp <- bp + geom_point() + geom_errorbarh(aes(xmax=max, xmin=min, height=.2))
	}else{	
		## order by 'Predictor', 'Val'
		df <- df[with(df,order(Predictor,-Val)),]
		df$Predictor <- factor(df$Predictor, levels=unique(df$Predictor))
		df$Method <- factor(df$Method, levels=rev(unique(df$Method)))
		bp <- ggplot(df, aes(Val, Method, colour=Predictor))
		bp <- bp + geom_point()
	}
	
	bp <- bp + scale_color_manual(values=xColormap("ggplot2")(length(levels(df$Predictor))))
	bp <- bp  + theme_bw() + theme(legend.position="right", legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=10,color="black"), axis.title.x=element_text(size=14,color="black"))
	bp <- bp + xlab(xlab)
	
	## caption
    if(signature){
    	caption <- paste("Created by xMLplot from Pi version", utils ::packageVersion("Pi"))
    	bp <- bp + labs(caption=caption) + theme(plot.caption=element_text(hjust=1,face='bold.italic',size=8,colour='#002147'))
    }
	
	## put arrows on x-axis
	bp <- bp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
	
	## x-axis position
	bp <- bp + scale_x_continuous(position="top")
	
	
	invisible(bp)
}
