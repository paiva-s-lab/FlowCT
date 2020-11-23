#' Transform a \code{fcs.SCE} object into a \code{flowset} object
#'
#' It tranforms a \code{fcs.SCE} object into a \href{https://www.rdocumentation.org/packages/flowCore/versions/1.38.2/topics/flowSet-class}{\code{flowset}} object.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which obtain a \code{flowset} object.
#' @keywords fcs.SCE to flowset
#' @keywords flowset generation
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @importFrom premessa as_flowFrame
#' @importFrom methods as
#' @examples

as.flowSet.SE <- function(fcs.SCE, assay.i){
  df_flowset <- lapply(unique(fcs.SCE$filename), function(x) as_flowFrame(t(assay(fcs.SCE[,fcs.SCE$filename == x], i = assay.i))))
  names(df_flowset) <- fcs.SCE@metadata$input_fcs
  df_flowset <- as(df_flowset, "flowSet")
  return(df_flowset)
}
