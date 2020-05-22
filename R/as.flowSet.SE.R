#'as.flowSet.SE
#'
#' It tranforms a \code{fcs.SCE} object into a \href{https://www.rdocumentation.org/packages/flowCore/versions/1.38.2/topics/flowSet-class}{\code{flowset}} object.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which obtain a \code{flowset} object.
#' @keywords fcs.SCE to flowset
#' @keywords flowset generation
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @importFrom premessa as_flowFrame
#' @importFrom methods as
#' @examples

as.flowSet.SE <- function(fcs.SCE, assay.i){
  df_flowset <- list()
  for(i in unique(fcs.SCE$filename)){
    aux_fcsSE <- fcs.SCE[,fcs.SCE$filename == i]
    df_flowset[[i]] <- as_flowFrame(t(assay(aux_fcsSE, i = assay.i)))
  }
  names(df_flowset) <- fcs.SCE@metadata$input_fcs
  df_flowset <- as(df_flowset, "flowSet")
  return(df_flowset)
}
