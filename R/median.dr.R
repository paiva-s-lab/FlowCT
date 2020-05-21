# 'median.dr
#'
#' This function draws a plot from internally calculated dimensional reduction (DR).
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param dr.method DR method to calculate. Possible values are "PCA", "tSNE" or "UMAP", individually, or "all". Take into account the number of samples for this DR, for few samples t-SNE and UMAP are useless.
#' @param perplexity.tsne Value for perplexity parameter in tSNE calculation (\href{https://distill.pub/2016/misread-tsne/}{more information}). For little number of samples this DR (t-SNE) cannot be calculated because preplexity restrictions. Default = \code{15}.
#' @param n.neighbors.umap Value for neighbors parameter in UMAP calculation (\href{https://www.math.upenn.edu/~jhansen/2018/05/04/UMAP/}{more information}). Default = \code{50}.
#' @param color.by Variable name (from \code{colData(fcs.SCE)}) for lines coloring. 
#' @param label.by Variable name (from \code{colData(fcs.SCE)}) for labelling. 
#' @param size Point size. Default = \code{3}.
#' @keywords dimensional reduction
#' @keywords median expression values
#' @export median.dr
#' @importFrom SummarizedExperiment assay
#' @importFrom utils capture.output
#' @examples
#' \dontrun{
#' median.dr(fcs, color.by = "condition")
#' }

median.dr <- function(fcs.SCE, assay.i = "normalized", markers.to.use = "all", dr.method = "PCA", perplexity.tsne = 15, n.neighbors.umap = 50, color.by, label.by = NULL, size = 3){
  data <- t(assay(fcs.SCE, i = assay.i))
  metadata <- fcs.SCE@metadata$reduced_metadata
  
  if(markers.to.use == "all") markers.to.use <- colnames(data)
  
  ## prepare median tables
  med <- median.values(fcs.SCE, assay.i = assay.i)
  
  ## dr
  if(!is.null(dr.method)){
    if(grepl("PCA|tSNE|UMAP", dr.method)){
      invisible(capture.output(dr <- dim.reduction(med, dr.method = dr.method, markers.to.use = markers.to.use, 
                          perplexity.tsne = perplexity.tsne, n.neighbors.umap = n.neighbors.umap)))
      
      ## plotting
      drmd <- merge(metadata, dr[[1]], by = "row.names")[-1]
      dr.plotting(drmd, plot.dr = dr.method, color.by = color.by, label.by = label.by, 
                        size = size, raster = F)
    }else{
      stop("Please, specify a valid reduction method: PCA, tSNE or UMAP.", call. = F)
    }
  }
}
