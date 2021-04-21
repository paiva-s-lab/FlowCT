#' Combine multiple subclustering with initial one
#'
#' It combines initial \code{fcs.SCE} object (without subclustering) with other \code{fcs.SCE} objects with subclustering analysis coming from downstream steps and generates a new \code{fcs.SCE} object. This final \code{fcs.SCE} object has an additional column combining all information from initial and subclustering analysis.
#' @param global.fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}. The initial one, without extracting any cell population.
#' @param subclustering.fcs.SCE A list with all \code{fcs.SCE} object(s) generated in the subclustering analysis (they have to come from the original \code{initial.fcs.SCE}).
#' @param clusters.named Column names from the \code{global.fcs.SCE} and \code{global.fcs.SCE} objects which contains renamed clusters (through \code{\link[FlowCT:clusters.rename]{FlowCT::clusters.rename()}}).
#' @keywords final fcs.SCE object
#' @keywords combine subclustering
#' @export
#' @importFrom SingleCellExperiment colData
#' @import dplyr
#' @examples
#' \dontrun{
#'  fcs_final <- combine.subclusterings(global.fcs.SCE = fcs, 
#'    clusters.named = c("SOM_named", "SOM_named_lymhpos", "SOM_named_monos"), 
#'    subclustering.fcs.SCE = list(fcs_lymphos, fcs_monos))
#' }

combine.subclusterings <- function(global.fcs.SCE, subclustering.fcs.SCE, clusters.named){
  mdg <- colData(global.fcs.SCE)
 
  subclusterings <- c()
  rm.cells <- list()
	  for(i in subclustering.fcs.SCE){
	    sub_pop <- as.character(unique(i[[clusters.named[1]]]))
	   
	    md_sub <- colData(i)
	   
	    # add differential cols from subclusterings
	    diff1 <- setdiff(colnames(md_sub), colnames(mdg))
	    mdg[,diff1] <- NA
	    diff2 <- setdiff(colnames(mdg), colnames(md_sub))
	    md_sub[,diff2] <- NA
	   
	    # extract those subclustered samples from original fcs.SCE object
	    mdg <- rbind(mdg[mdg[[clusters.named[1]]] != sub_pop,], md_sub)

	    # delete removed populations from the (i)th subclustering within the first fcs.SCE object
	    if(!is.null(i@metadata$removed_populations)){
	      for(j in names(i@metadata$removed_populations)){
	          rm.cells[[j]] <- i@metadata$removed_populations[[j]]
	          mdg <- mdg[setdiff(mdg$cell_id, i@metadata$removed_populations[[j]]),]
	         
	          # substract deleted cells from original assays
	          global.fcs.SCE <- global.fcs.SCE[, global.fcs.SCE$cell_id %in% setdiff(mdg$cell_id, i@metadata$removed_populations[[j]])]
	       }
	    }
	  }
 
  # create final named_cluster col
  for(i in setdiff(colnames(mdg), c(colnames(colData(global.fcs.SCE)), clusters.named)))
  	mdg[,i] <- ifelse(is.na(mdg[,i]), 0, mdg[,i]) #replace NAs by 0 to avoid FCS wrong building
 
  mdg$final_clustering <- factor(do.call(dplyr::coalesce, mdg[,clusters.named[-1]]))
  mdg$final_clustering <- ifelse(is.na(mdg$final_clustering), as.character(mdg[[clusters.named[1]]]), as.character(mdg$final_clustering))
 
  colData(global.fcs.SCE) <- mdg[match(global.fcs.SCE$cell_id, mdg$cell_id),] # same cell_id order than initial fcs.SCE
  global.fcs.SCE@metadata$subclusterings$populations <- paste(subclusterings, collapse = " + ")
  global.fcs.SCE@metadata$subclusterings$removed_populations <- rm.cells
 
  return(global.fcs.SCE)
}
