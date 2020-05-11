# 'combine.subclusterings
#'
#' It combines initial \code{FCS.SE} object (without subclustering) with other \code{FCS.SE} objects with subclustering analysis coming from downstream steps and generates a new \code{FCS.SE} object. This final \code{FCS.SE} object has an additional column combining all information from initial and subclustering analysis.
#' @param initial.fcs.SE A \code{FCS.SE} object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}. The initial one, without extracting any cell population.
#' @param subclustering.fcs.SE A list with all \code{FCS.SE} object(s) generated in the subclustering analysis (they have to come from the original \code{initial.fcs.SE}).
#' @param clusters.named Column name from the \code{initial.fcs.SE} object which contains renamed clusters (through \code{\link[FlowCT:clusters.rename]{FlowCT::clusters.rename()}}) and has been used to extract cell populations for subclustering steps.
#' @keywords final FCS.SE object
#' @keywords combine subclustering
#' @export
#' @importFrom SummarizedExperiment colData
#' @examples
#' \dontrun{
#'  fcs_se_final <- combine.subclusterings(initial.fcs.SE = fcs_se, clusters.named = "SOM_named", subclustering.fcs.SE = list(fcs_se_lymphos, fcs_se_monos))
#' }

combine.subclusterings <- function(initial.fcs.SE, subclustering.fcs.SE, clusters.named = "SOM_named"){
  mdg <- colData(initial.fcs.SE)
  
  subclusterings <- c()
  rm.cells <- list()
  for(i in subclustering.fcs.SE){
    subclusterings <- c(subclusterings, i@metadata$subclustering)
    
    md_sub <- colData(i)
    
    # add differential cols from subclusterings
    diff1 <- setdiff(colnames(md_sub), colnames(mdg))
    mdg[,diff1] <- NA
    diff2 <- setdiff(colnames(mdg), colnames(md_sub))
    md_sub[,diff2] <- NA
    
    # extract those subclustered samples from original fcs.SE object
    mdg <- rbind(mdg[mdg[,clusters.named] != i@metadata$subclustering,], md_sub)
    
    # delete removed populations from the (i)th subclustering within the first fcs.SE object
    if(!is.null(i@metadata$removed_populations)){
      for(j in names(i@metadata$removed_populations))
        rm.cells[[j]] <- i@metadata$removed_populations[[j]]
        mdg <- mdg[setdiff(mdg$cell_id, i@metadata$removed_populations[[j]]),]
        
        # substract deleted cells from original assays
        initial.fcs.SE <- initial.fcs.SE[, initial.fcs.SE$cell_id %in% setdiff(mdg$cell_id, i@metadata$removed_populations[[j]])]
    }
  }
  
  # create final named_cluster col ---> beta, this will bring problems with multiple subclusterings!
  diff <- setdiff(colnames(mdg), colnames(colData(initial.fcs.SE)))
  for(i in diff){
    if(suppressWarnings(sum(is.na(as.numeric(as.character(mdg[,i])))) == nrow(mdg))){ #detect names_clusters cols in subclustering
      mdg$tmp <- ifelse(is.na(mdg[,i]), as.character(mdg[,clusters.named]), 
                        as.character(mdg[,i]))
    }
    mdg[,i] <- ifelse(is.na(mdg[,i]), 0, mdg[,i]) #replace NAs by 0 to avoid FCS wrong building
  }
  colnames(mdg)[ncol(mdg)] <- paste0(clusters.named, "_final")
  
  SummarizedExperiment::colData(initial.fcs.SE) <- mdg
  initial.fcs.SE@metadata$subclusterings$populations <- paste(subclusterings, collapse = " + ")
  initial.fcs.SE@metadata$subclusterings$removed_populations <- rm.cells
  
  return(initial.fcs.SE)
}