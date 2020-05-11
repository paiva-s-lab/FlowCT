#'as.flowSet.SE
#'
#' It tranforms a FCS.SE object into a (\href{https://www.rdocumentation.org/packages/flowCore/versions/1.38.2/topics/flowSet-class}{flowset}) object.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which obtain a \code{flowset} object.
#' @keywords FCS.SE to flowset
#' @keywords flowset generation
#' @export
#' @importFrom SummarizedExperiment colData
#' @examples
#' \dontrun{
#'  sc_metadata <- scMetadata.fromFCS(fcs, metadata, add.exprs = T)
#' }

as.flowSet.SE <- function(fcs.SE, assay.i){
  df_flowset <- list()
  for(i in unique(fcs.SE$filename)){
    aux_fcsSE <- fcs.SE[,fcs.SE$filename == i]
    df_flowset[[i]] <- premessa::as_flowFrame(t(assay(aux_fcsSE, i = "raw")))
  }
  names(df_flowset) <- fcs.SE@metadata$input_fcs
  df_flowset <- methods::as(df_flowset, "flowSet")
  return(df_flowset)
}
