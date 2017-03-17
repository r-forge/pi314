#' Function to visualise cross-validation performance against tuning parameters
#'
#' \code{xMLcompare} is supposed to visualise cross-validation performance against tuning parameters. 
#'
#' @param list_ML a list of class "train" or "train.formula" objects (resulting from caret::train)
#' @param metric the performance metric to plot. It can be one of 'ROC', 'Sens' (Sensitivity) and 'Spec' (Specificity)
#' @param xlab a title for the x axis
#' @param xlimits the limit for the x axis
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{xMLcompare}}
#' @include xMLcompare.r
#' @examples
#' \dontrun{
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"
#' \dontrun{
# ls_ML <- list(GBM=fit_gbm, SVM=fit_svm, RDA=fit_rda, KNN=fit_knn, PLS=fit_pls, NNET=fit_nnet, RF=fit_myrf, CRF=fit_crf, GLM=fit_glm, BGLM=fit_bglm, BLR=fit_blr)
# names(ls_ML) <- c("Gradient Boosting Machine (GBM)", "Support Vector Machine (SVM)", "Regularized Discriminant Analysis (RDA)", "K-Nearest Neighbor (KNN)", "Partial Least Squares (PLS)", "Neural Network (NNET)", "Random Forest (RF)", "Conditional Inference Random Forest (CIRF)", "Generalized Linear Model (GLM)", "Bayes Generalized Linear Model (BGLM)", "Boosted Logistic Regression (BLR)")
#' gp <- xMLcompare(ls_ML, xlimits=c(0.5,1))
#' }

xMLcompare <-function(list_ML, metric=c("ROC","Sens","Spec"), xlab=NA, xlimits=c(0.5,1))
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    metric <- match.arg(metric)

    if(all(class(list_ML) %in% c("train","train.formula"))){
    	list_ML <- list(list_ML)
    }else if(class(list_ML)=="list"){
		## Remove null elements in a list
		list_ML <- base::Filter(base::Negate(is.null), list_ML)
		if(length(list_ML)==0){
			warnings("All train objects are NULL!")
			return(NULL)
		}
    }else{
    	stop("The function must apply to a 'list' object, or a 'train'/'train.formula' object.\n")
    }
	
	resamps <- caret::resamples(list_ML)
	if(metric=="ROC"){
		ind <- grep("ROC",colnames(resamps$values))
		df_ML <- resamps$values[,ind]
		colnames(df_ML) <- gsub("~ROC", "", colnames(df_ML))
		if(is.na(xlab)){
			xlab <- "AUC: 95% Confidence Interval\n(Repeated Cross-Validation)"
		}
		
	}else if(metric=="Sens"){
		ind <- grep("Sens",colnames(resamps$values))
		df_ML <- resamps$values[,ind]
		colnames(df_ML) <- gsub("~Sens", "", colnames(df_ML))
		if(is.na(xlab)){
			xlab <- "Sensitivity: 95% Confidence Interval\n(Repeated Cross-Validation)"
		}
		
	}else if(metric=="Spec"){
		ind <- grep("Spec",colnames(resamps$values))
		df_ML <- resamps$values[,ind]
		colnames(df_ML) <- gsub("~Spec", "", colnames(df_ML))
		if(is.na(xlab)){
			xlab <- "Specificity: 95% Confidence Interval\n(Repeated Cross-Validation)"
		}
		
	}
	
	res_ttest <- lapply(df_ML, t.test)
	ls_res <- lapply(1:length(res_ttest), function(i){
		x <- res_ttest[[i]]
		data.frame(ML=names(res_ttest)[i], mean=x$estimate, conf_lower=x$conf.int[1], conf_upper=x$conf.int[2])
	})
	df <- do.call(rbind, ls_res)
	rownames(df) <- NULL
	
	## order by 'mean'
	mean <- ML <- ''
	df <- df[with(df,order(-mean)),]
	df$ML <- factor(df$ML, levels=rev(unique(df$ML)))
	
	conf_upper <- conf_lower <- ''
	bp <- ggplot(df, aes(mean, ML))
	bp <- bp + geom_point() + geom_errorbarh(aes(xmax=conf_upper, xmin=conf_lower, height=.2))
	bp <- bp + scale_color_manual(values=xColormap("ggplot2")(length(levels(df$method))))
	bp <- bp  + theme_bw() + theme(legend.position="right", legend.title=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=10,color="black"), axis.title.x=element_text(size=14,color="black"), panel.background=element_rect(fill=rgb(0.95,0.95,0.95,1)))
	bp <- bp + xlab(xlab)
	## put arrows on x-axis
	bp <- bp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
	gp <- bp + scale_x_continuous(position="top", limits=xlimits)
	
	## append 'CI' and 'Para'
	### CI
	df_gp <- gp$data
	df_gp$CI <- paste0('[', signif(df_gp$conf_lower,digits=3), ', ', signif(df_gp$conf_upper,digits=3), ']')
	### Para
	ind <- match(df_gp$ML, names(list_ML))
	ls_ML <- list_ML[ind]
	df_gp$Para <- sapply(ls_ML, function(x){
		tmp <- x$bestTune
		if(all(tmp!='none')){
			paste0(names(tmp),'=',tmp, collapse='; ')
		}else{
			'None'
		}
	})
	gp$data <- df_gp
	
    invisible(gp)
}

    
