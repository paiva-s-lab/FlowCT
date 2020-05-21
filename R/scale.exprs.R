#' scale.exprs
#'
#' It scales a matrix expression and fit all values to a range from 0 to 1.
#' @param data A fcs.SCE object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}} or a expression table with events in rows and markers in columns.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param include.FCS.SCE Logical indicating if new scaled matrix must to be stored within the \code{fcs.SCE} object.
#' @param scaled.matrix.name New scaled matrix name (it will stored within the \code{fcs.SCE} object). Default = \code{"scaled"}.
#' @keywords scale
#' @export
#' @importFrom matrixStats colQuantiles
#' @importFrom SummarizedExperiment assay
#' @examples
#' \dontrun{
#' # option 1: save scaled matrix within the fcs.SCE object (if provided)
#' fcs <- scale.exprs(fcs.SCE = fcs)
#'
#' # option 2: scale an external data.frame
#' markers_expresion_scaled <- scale.exprs(fcs.SCE = markers_expresion)
#' }

scale.exprs <- function(data, assay.i = "normalized", include.FCS.SCE = F, scaled.matrix.name = "scaled"){
	if(class(data) == "SingleCellExperiment") data <- t(assay(data, i = assay.i))
	rng <- colQuantiles(data, probs = c(0.01, 0.99))
	expr01 <- t((t(data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
	expr01[expr01 < 0] <- 0
	expr01[expr01 > 1] <- 1


	if(class(data) == "SingleCellExperiment" & include.FCS.SCE){
		assay(fcs.SE, i = scaled.matrix.name) <- expr01
		return(fcs.SE)
	}else{
		return(expr01)
	}

}
