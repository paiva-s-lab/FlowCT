# 'median.values
#'
#' It calculates median values according a specfied variable, normally "filename".
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate medians. Default = \code{"raw"}.
#' @param var Variable for grouping and calculating medians. Default = \code{"filename"} (i.e., names of FCS files).
#' @keywords median values
#' @keywords MFI
#' @keywords median fluorescence intensity
#' @export median.values
#' @examples
#' \dontrun{
#' mfis_FCS <- median.values(fcs.SE = fcs_se)
#' med_SOM_clust <- median.values(fcs.SE = fcs_se, assay.i = "normalized", var = "SOM_named")
#' }

median.values <- function(fcs.SE, assay.i = "raw", var = "filename"){
  med <- list()
  for(i in unique(fcs.SE[[var]])){
    aux_se <- fcs.SE[,fcs.SE[[var]] == i]
    med[[i]] <- apply(SummarizedExperiment::assay(aux_se, i = assay.i), 1, stats::median)
  }
  if(is.null(names(med))) names(med) <- unique(fcs.SE[[var]])
  return(t(as.data.frame(med)))
}
