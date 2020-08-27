#' Normalization for markers' expression
#'
#' It normalizes expression abnormalities within selected markers.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"transformed"}.
#' @param marker Marker name(s) to normalize (based on visual inspection with \code{\link[FlowCT.v2:multidensity]{FlowCT.v2::multidensity()}}).
#' @param method Methodology to perform each marker normalization. Possible values are "gauss" or "warp" (from \code{flowStats} package).
#' @param new.matrix.name New normalized matrix name (it will stored within the \code{fcs.SCE} object). Default = \code{"normalized"}.
#' @keywords marker alignment
#' @keywords marker normalization
#' @keywords gaussian warp
#' @export marker.normalization
#' @importFrom flowStats gaussNorm warpSet
#' @importFrom SummarizedExperiment assay
#' @importFrom flowCore fsApply exprs
#' @importFrom utils capture.output
#' @examples
#' \dontrun{
#'  fcs_gauss <- normalization.flw(fcs.SCE = fcs, marker = c("CD62L", "CCR4", "SSC_A"), method = "gauss")
#'  fcs_warp <- normalization.flw(fcs.SCE = fcs, marker = c("CD62L", "CCR4", "SSC_A"), method = "warp")
#' }

marker.normalization <- function(fcs.SCE, assay.i = "transformed", method, marker, new.matrix.name = "normalized"){
  fcs1 <- as.flowSet.SE(fcs.SCE, assay.i = assay.i)
  if(method == "gauss"){
    for(i in marker){
      cat("Applying normalization to:", i, "\n")
      tryCatch(
        {invisible(capture.output(fcs1 <- gaussNorm(fcs1, i, peak.density.thr = 0.001)$flowset))},
        error = function(e){
          message(paste0(" !Normalization adjusted to [1] landmark for this marker (", i, ")"))
          invisible(capture.output(fcs1 <- gaussNorm(fcs1, i, peak.density.thr = 0.001, max.lms = 1)$flowset))
        }
      )
    }
  }else if(method == "warp"){
    for(i in marker){
      cat("Applying normalization to:", i, "\n")
      invisible(capture.output(fcs1 <- warpSet(fcs1, i)))
    }
  }else{
    stop("Please, indicate a valid normalization method: gauss or warp.", call. = F)
  }
  assay(fcs.SCE, i = new.matrix.name) <- t(fsApply(fcs1, exprs))
  return(fcs.SCE)
}
