# 'remove.pop
#'
#' It removes one or multiple cell populations from a \code{fcs.SCE} object.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT]{fcs.SCE()}}.
#' @param population Name(s) of cell population(s) to be removed.
#' @param clusters.named Column name from the \code{colData(fcs.SCE)} object which contains renamed clusters (through \code{\link[FlowCT]{clusters.rename()}}).
#' @keywords remove cell population
#' @keywords debris
#' @export
#' @importFrom SummarizedExperiment colData
#' @examples
#' \dontrun{
#' fcs_rm <- remove.pop(fcs, clusters.named = "SOM_named", 
#'     population = c("debris", "unclassified"))
#' }

remove.pop <- function(fcs.SCE, population, clusters.named){
  for(i in population){
    fcs.SCE@metadata$removed_populations[[i]] <- colData(fcs.SCE)[fcs.SCE[[clusters.named]] == i,"cell_id"]
    fcs.SCE <- fcs.SCE[,fcs.SCE[[clusters.named]] != i]
    fcs.SCE[[clusters.named]] <- droplevels(fcs.SCE[[clusters.named]])
  }
  return(fcs.SCE)
}
