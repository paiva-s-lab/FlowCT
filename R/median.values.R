# 'median.values
#'
#' It calculates median values according a specfied variable, normaly "filename".
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate medians. Default = \code{"raw"}.
#' @param var Variable for grouping and calculating medians. Default = \code{"filename"} (i.e., names of FCS files).
#' @keywords median values
#' @keywords MFI
#' @keywords median fluorescence intensity
#' @export median.values
#' @importFrom SummarizedExperiment	assay
#' @importFrom stats median
#' @examples
#' \dontrun{
#' mfis_FCS <- median.values(fcs.SCE = fcs)
#' med_SOM_clust <- median.values(fcs.SCE = fcs, assay.i = "normalized", var = "SOM_named")
#' }

median.values <- function(fcs.SCE, assay.i = "raw", var = "filename"){
  med <- list()
  for(i in unique(fcs.SCE[[var]])){
    aux_se <- fcs.SCE[,fcs.SCE[[var]] == i]
    med[[i]] <- apply(assay(aux_se, i = assay.i), 1, median)
  }
  if(is.null(names(med))) names(med) <- unique(fcs.SCE[[var]])
  return(t(as.data.frame(med)))
}
