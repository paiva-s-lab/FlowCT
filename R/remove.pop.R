# 'remove.pop
#'
#' It removes one or multiple cell populations from a \code{FCS.SE} object.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param population Name(s) of cell population(s) to be removed.
#' @param clusters.named Column name from the \code{colData(fcs.SE)} object which contains renamed clusters (through \code{\link[FlowCT:clusters.rename]{FlowCT::clusters.rename()}}).
#' @keywords remove cell population
#' @keywords debris
#' @export
#' @importFrom SummarizedExperiment colData
#' @examples
#' \dontrun{
#' fcs_se_rm <- remove.pop(fcs_se, clusters.named = "SOM_named", 
#'     population = c("debris", "unclassified"))
#' }

remove.pop <- function(fcs.SE, population, clusters.named){
  for(i in population){
    fcs.SE@metadata$removed_populations[[i]] <- colData(fcs.SE)[fcs.SE[[clusters.named]] == i,"cell_id"]
    fcs.SE <- fcs.SE[,fcs.SE[[clusters.named]] != i]
    fcs.SE[[clusters.named]] <- droplevels(fcs.SE[[clusters.named]])
  }
  return(fcs.SE)
}
