#' Normalization for markers' expression
#'
#' It normalizes expression abnormalities within selected markers or within the entire expression set.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"transformed"}.
#' @param marker.to.norm Marker name(s) to normalize (based on visual inspection with \code{\link[FlowCT.v2:multidensity]{FlowCT.v2::multidensity()}}). Only applicable for methods "gauss" and "warp".
#' @param norm.method Methodology to perform each marker normalization. Possible values are "gauss" or "warp" (from \code{flowStats} package) or "harmony" (from \code{harmony} package for single-cell).
#' @param var.to.use Variable from \code{colData(fcs.SCE)} for calculating the normalization (only if \code{norm.method = "harmony}).
#' @param norm.matrix.name New normalized matrix name (it will stored within the \code{fcs.SCE} object). Default = \code{"normalized"}.
#' @keywords marker alignment
#' @keywords marker normalization
#' @keywords gaussian warp harmony
#' @export normalization.flw
#' @importFrom flowStats gaussNorm warpSet
#' @importFrom SummarizedExperiment assay
#' @importFrom flowCore fsApply exprs
#' @importFrom utils capture.output
#' @importFrom harmony HarmonyMatrix
#' @examples
#' \dontrun{
#'  fcs_gauss <- normalization.flw(fcs.SCE = fcs, marker.to.norm = c("CD62L", "CCR4", "SSC_A"), norm.method = "gauss")
#'  fcs_warp <- normalization.flw(fcs.SCE = fcs, marker.to.norm = c("CD62L", "CCR4", "SSC_A"), norm.method = "warp")
#'  fcs_harmony <- normalization.flw(fcs.SCE = fcs, norm.method = "harmony", var.to.use = "condition")
#' }

normalization.flw <- function(fcs.SCE, assay.i = "transformed", marker.to.norm, norm.method, var.to.use, norm.matrix.name = "normalized"){
  fcs <- as.flowSet.SE(fcs.SCE, assay.i = assay.i)
  if(norm.method == "gauss"){
    for(i in marker.to.norm){
      cat("Applying normalization to:", i, "\n")
      tryCatch(
        {invisible(capture.output(data <- gaussNorm(fcs, i, peak.density.thr = 0.001)))},
        error = function(e){
          message(paste0(" !Normalization adjusted to [1] landmark for this marker (", i, ")"))
          invisible(capture.output(data <- gaussNorm(fcs, i, peak.density.thr = 0.001, max.lms = 1)))
        }
      )
    }
    norm_data <- t(fsApply(data, exprs))
  }else if(norm.method == "warp"){
    for(i in marker.to.norm){
      cat("Applying normalization to:", i, "\n")
      invisible(capture.output(data <- warpSet(fcs, i)))
    }
    norm_data <- t(fsApply(data, exprs))
  }else if(norm.method == "harmony"){
    norm_data <- HarmonyMatrix(data_mat = assay(fcs.SCE, assay.i), meta_data = colData(fcs.SCE),
                               vars_use = var.to.use, do_pca = F, verbose = F)
  }else{
    cat("Please, indicate a valid normalization method: 'gauss' or 'warp'.")
  }
  assay(fcs.SCE, i = norm.matrix.name) <- norm_data
  return(fcs.SCE)
}
