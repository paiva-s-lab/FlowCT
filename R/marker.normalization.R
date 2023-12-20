#' Normalization for markers' expression
#'
#' It normalizes expression abnormalities within selected markers.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"transformed"}.
#' @param marker Marker name(s) to normalize (based on visual inspection with \code{\link[FlowCT:multidensity]{FlowCT::multidensity()}}).
#' @param method Methodology to perform each marker normalization. Possible values are "gauss" or "warp" (from \code{flowStats} package).
#' @param new.matrix.name New normalized matrix name (it will stored within the \code{fcs.SCE} object). Default = \code{"normalized"}.
#' @keywords marker alignment
#' @keywords marker normalization
#' @keywords gaussian warp
#' @export marker.normalization
#' @importFrom SummarizedExperiment assay
#' @importFrom flowCore fsApply exprs
#' @importFrom utils capture.output
#' @examples
#' \dontrun{
#'  fcs_gauss <- normalization.flw(fcs.SCE = fcs, method = "gauss", 
#'      marker = c("CD62L", "CCR4", "SSC_A"))
#'  fcs_warp <- normalization.flw(fcs.SCE = fcs, method = "warp", 
#'      marker = c("CD62L", "CCR4", "SSC_A"))
#' }

marker.normalization <- function(fcs.SCE, assay.i = "transformed", method, marker, new.matrix.name = "normalized"){
  if (!requireNamespace("flowStats", quietly = TRUE)) cat("Package \"flowStats\" is needed for this function to work. Installing it...\n\n")
  BiocManager::install("flowStats")
  
  aux_fcs <- as.flowSet.SE(fcs.SCE, assay.i = assay.i)
  if(method == "gauss"){
    for(i in marker){
      cat("Applying normalization to:", i, "\n")
      tryCatch(
        {invisible(capture.output(aux_fcs <- flowStats::gaussNorm(aux_fcs, i, peak.density.thr = 0.001)$flowset))},
        error = function(e){
          message(paste0(" !Normalization adjusted to [1] landmark for this marker (", i, ")"))
          invisible(capture.output(aux_fcs <- flowStats::gaussNorm(aux_fcs, i, peak.density.thr = 0.001, max.lms = 1)$flowset))
        }
      )
    }
  }else if(method == "warp"){
    for(i in marker){
      cat("Applying normalization to:", i, "\n")
      invisible(capture.output(aux_fcs <- flowStats::warpSet(aux_fcs, i)))
    }
  }else{
    stop("Please, indicate a valid normalization method: gauss or warp.", call. = F)
  }
  
  norm_matrix <- t(flowCore::fsApply(aux_fcs, exprs))

  assay(fcs.SCE, i = new.matrix.name, withDimnames = F) <- norm_matrix
  return(fcs.SCE)
}
