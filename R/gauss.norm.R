#'gauss.norm
#'
#' It performs a Gaussian normalization based on \code{\link[flowStats:gaussNorm]{flowStats::gaussNorm()}} function.
#' @param fcs.SCE A fcs.SCE object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"transformed"}.
#' @param marker.to.norm Marker name(s) to normalize (based on visual inspection with \code{\link[FlowCT:multidensity]{FlowCT::multidensity()}}).
#' @param norm.matrix.name New normalized matrix name (it will stored within the \code{fcs.SCE} object). Default = \code{"normalized"}.
#' @keywords Gaussian normalization
#' @keywords marker alignment
#' @keywords marker normalization
#' @export
#' @importFrom flowStats gaussNorm
#' @importFrom SummarizedExperiment assay
#' @importFrom flowCore fsApply exprs
#' @examples
#' \dontrun{
#'  fcs <- gauss.norm(fcs.SCE = fcs, marker.to.norm = c("CD62L", "CCR4", "SSC_A"))
#' }

gauss.norm <- function(fcs.SCE, assay.i = "transformed", marker.to.norm, norm.matrix.name = "normalized"){
  data <- as.flowSet.SE(fcs.SCE, assay.i = assay.i)
  for(i in marker.to.norm){
    data <- gaussNorm(data, i)$flowset
  }
  assay(fcs.SCE, i = norm.matrix.name) <- t(fsApply(data, exprs))
  return(fcs.SCE)
}
