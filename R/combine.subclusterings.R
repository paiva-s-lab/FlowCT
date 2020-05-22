# 'combine.subclusterings
#'
#' It combines initial \code{fcs.SCE} object (without subclustering) with other \code{fcs.SCE} objects with subclustering analysis coming from downstream steps and generates a new \code{fcs.SCE} object. This final \code{fcs.SCE} object has an additional column combining all information from initial and subclustering analysis.
#' @param initial.fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{fcs.SCE()}}. The initial one, without extracting any cell population.
#' @param subclustering.fcs.SCE A list with all \code{fcs.SCE} object(s) generated in the subclustering analysis (they have to come from the original \code{initial.fcs.SCE}).
#' @param clusters.named Column name from the \code{initial.fcs.SCE} object which contains renamed clusters (through \code{\link[FlowCT:clusters.rename]{clusters.rename()}}) and has been used to extract cell populations for subclustering steps.
#' @keywords final fcs.SCE object
#' @keywords combine subclustering
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData
#' @examples
#' \dontrun{
#'  fcs_final <- combine.subclusterings(initial.fcs.SCE = fcs, clusters.named = "SOM_named", subclustering.fcs.SCE = list(fcs_lymphos, fcs_monos))
#' }

combine.subclusterings <- function(initial.fcs.SCE, subclustering.fcs.SCE, clusters.named = "SOM_named"){
  mdg <- colData(initial.fcs.SCE)
  
  subclusterings <- c()
  rm.cells <- list()
  for(i in subclustering.fcs.SCE){
    subclusterings <- c(subclusterings, i@metadata$subclustering)
    
    md_sub <- colData(i)
    
    # add differential cols from subclusterings
    diff1 <- setdiff(colnames(md_sub), colnames(mdg))
    mdg[,diff1] <- NA
    diff2 <- setdiff(colnames(mdg), colnames(md_sub))
    md_sub[,diff2] <- NA
    
    # extract those subclustered samples from original fcs.SCE object
    mdg <- rbind(mdg[mdg[,clusters.named] != i@metadata$subclustering,], md_sub)
    
    # delete removed populations from the (i)th subclustering within the first fcs.SCE object
    if(!is.null(i@metadata$removed_populations)){
      for(j in names(i@metadata$removed_populations))
        rm.cells[[j]] <- i@metadata$removed_populations[[j]]
        mdg <- mdg[setdiff(mdg$cell_id, i@metadata$removed_populations[[j]]),]
        
        # substract deleted cells from original assays
        initial.fcs.SCE <- initial.fcs.SCE[, initial.fcs.SCE$cell_id %in% setdiff(mdg$cell_id, i@metadata$removed_populations[[j]])]
    }
  }
  
  # create final named_cluster col ---> beta, this will bring problems with multiple subclusterings!
  diff <- setdiff(colnames(mdg), colnames(colData(initial.fcs.SCE)))
  for(i in diff){
    if(suppressWarnings(sum(is.na(as.numeric(as.character(mdg[,i])))) == nrow(mdg))){ #detect names_clusters cols in subclustering
      mdg$tmp <- ifelse(is.na(mdg[,i]), as.character(mdg[,clusters.named]), 
                        as.character(mdg[,i]))
    }
    mdg[,i] <- ifelse(is.na(mdg[,i]), 0, mdg[,i]) #replace NAs by 0 to avoid FCS wrong building
  }
  colnames(mdg)[ncol(mdg)] <- paste0(clusters.named, "_final")
  
  colData(initial.fcs.SCE) <- mdg
  initial.fcs.SCE@metadata$subclusterings$populations <- paste(subclusterings, collapse = " + ")
  initial.fcs.SCE@metadata$subclusterings$removed_populations <- rm.cells
  
  return(initial.fcs.SCE)
}