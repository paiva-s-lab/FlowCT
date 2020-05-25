#'normalization.flw
#'
#' It normalizes expression abnormalities within selected markers.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"transformed"}.
#' @param norm.method Methodology to perform each marker normalization. Possible values are \code{gauss} and \code{warp}.
#' @param marker.to.norm Marker name(s) to normalize (based on visual inspection with \code{\link[FlowCT.v2:multidensity]{FlowCT.v2::multidensity()}}).
#' @param norm.matrix.name New normalized matrix name (it will stored within the \code{fcs.SCE} object). Default = \code{"normalized"}.
#' @keywords marker alignment
#' @keywords marker normalization
#' @export
#' @importFrom flowStats gaussNorm warpSet
#' @importFrom SummarizedExperiment assay
#' @importFrom flowCore fsApply exprs
#' @importFrom utils capture.output
#' @examples
#' \dontrun{
#'  fcs_gauss <- normalization.flw(fcs.SCE = fcs, marker.to.norm = c("CD62L", "CCR4", "SSC_A"), norm.method = "gauss")
#'  fcs_warp <- normalization.flw(fcs.SCE = fcs, marker.to.norm = c("CD62L", "CCR4", "SSC_A"), norm.method = "warp")
#' }

normalization.flw <- function(fcs.SCE, assay.i = "transformed", marker.to.norm, norm.method, norm.matrix.name = "normalized"){
  fcs <- as.flowSet.SE(fcs.SCE, assay.i = assay.i)
  if(norm.method == "gauss"){
    for(i in marker.to.norm){
      cat("Applying normalization to:", i, "\n")
      tryCatch(
        {invisible(capture.output(data <- gaussNorm(fcs, i, peak.density.thr = 0.001)$flowset))},
        error = function(e){
          message(paste0(" !Normalization adjusted to [1] landmark for this marker (", i, ")"))
          invisible(capture.output(data <- gaussNorm(fcs, i, peak.density.thr = 0.001, max.lms = 1)$flowset))
          # return(data)
        }
      )
    }
  }else if(norm.method == "warp"){
    for(i in marker.to.norm){
      cat("Applying normalization to:", i, "\n")
      invisible(capture.output(data <- warpSet(fcs, i, clipRange = 0.5)$flowset))
    }
  }else{
    cat("Please, indicate a valid normalization method: 'gauss' or 'warp'.")
  }
  assay(fcs.SCE, i = norm.matrix.name) <- t(fsApply(data, exprs))
  return(fcs.SCE)
}
