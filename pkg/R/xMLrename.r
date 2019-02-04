#' Function to rename predictors used in machine learning
#'
#' \code{xMLrename} is supposed to rename predictors used in machine learning. It returns an object of class "sTarget".
#'
#' @param sTarget an object of class "sTarget"
#' @param old_names a vector for the original names of predictors to be renamed
#' @param new_names a vector for the new names
#' @return an object of class "sTarget"
#' @note none
#' @export
#' @seealso \code{\link{xMLrandomforest}}
#' @include xMLrename.r
#' @examples
#' \dontrun{
#' # Load the library
#' library(Pi)
#' }
#' RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
#' \dontrun{
#' old_names <- colnames(sTarget$predictor)[-c(1,2)]
#' new_names <- c('nearbyGenes_nGene: nearby genes', 'eQTL_eGene: resting state (CD14+)', 'eQTL_eGene: activating state (CD14+ by LPS2h)', 'eQTL_eGene: activating state (CD14+ by LPS24h)', 'eQTL_eGene: activating state (CD14+ by IFN24h)', 'eQTL_eGene: B cells', 'eQTL_eGene: CD4+ T cells', 'eQTL_eGene: CD8+ T cells', 'eQTL_eGene: neutrophils', 'eQTL_eGene: NK cells','eQTL_eGene: peripheral blood', 'HiC_cGene: monocytes', 'HiC_cGene: macrophages (M0)', 'HiC_cGene: macrophages (M1)', 'HiC_cGene: macrophages (M2)', 'HiC_cGene: neutrophils', 'HiC_cGene: CD4+ T cells (naive)', 'HiC_cGene: CD4+ T cells (total)', 'HiC_cGene: CD8+ T cells (naive)', 'HiC_cGene: CD8+ T cells (total)', 'HiC_cGene: B cells (naive)', 'HiC_cGene: B cells (total)', 'Annotation_dGene: disease genes', 'Annotation_pGene: phenotype genes', 'Annotation_fGene: function genes')
#' sTarget_rename <- xMLrename(sTarget, old_names, new_names)
#' }

xMLrename <- function(sTarget, old_names, new_names) 
{
    
    if(class(sTarget) != "sTarget"){
    	stop("The function must apply to a 'sTarget' object.\n")
    }
    sTarget_renamed <- sTarget
    
    ## predictor
	if(!is.null(sTarget_renamed$predictor)){
		ind <- match(colnames(sTarget_renamed$predictor), old_names)
		colnames(sTarget_renamed$predictor)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
	
	## importance2fold
	if(!is.null(sTarget_renamed$importance2fold)){
		ind <- match(rownames(sTarget_renamed$importance2fold), old_names)
		rownames(sTarget_renamed$importance2fold)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
		
	## roc2fold
	if(!is.null(sTarget_renamed$roc2fold)){
		ind <- match(rownames(sTarget_renamed$roc2fold), old_names)
		rownames(sTarget_renamed$roc2fold)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
	
	## fmax2fold
	if(!is.null(sTarget_renamed$fmax2fold)){
		ind <- match(rownames(sTarget_renamed$fmax2fold), old_names)
		rownames(sTarget_renamed$fmax2fold)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
	
	## importance
	if(!is.null(sTarget_renamed$importance)){
		ind <- match(rownames(sTarget_renamed$importance), old_names)
		rownames(sTarget_renamed$importance)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
	
	## performance
	if(!is.null(sTarget_renamed$performance)){	
		ind <- match(rownames(sTarget_renamed$performance), old_names)
		rownames(sTarget_renamed$performance)[!is.na(ind)] <- new_names[ind[!is.na(ind)]]
	}
	
	invisible(sTarget_renamed)
}
