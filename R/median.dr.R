#' Calculate a DR with median values
#'
#' This function draws a plot from internally calculated dimensional reduction (DR).
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param dr.method DR method to calculate. Possible values are "PCA", "tSNE", "UMAP" or "kmeans". Take into account the number of samples for this DR, for few samples t-SNE and UMAP are useless. Default = \code{"kmeans"}.
#' @param num.k Number of clusters if \code{dr.method = kmeans}. Default = \code{3}.
#' @param perplexity.tsne Value for perplexity parameter in tSNE calculation (\href{https://distill.pub/2016/misread-tsne/}{more information}). For little number of samples this DR (t-SNE) cannot be calculated because preplexity restrictions. Default = \code{15}.
#' @param n.neighbors.umap Value for neighbors parameter in UMAP calculation (\href{https://www.math.upenn.edu/~jhansen/2018/05/04/UMAP/}{more information}). Default = \code{50}.
#' @param color.by Variable name (from \code{colData(fcs.SCE)}) for lines coloring. Default = \code{"filename"}.
#' @param label.by Variable name (from \code{colData(fcs.SCE)}) for labelling. Default = \code{NULL}.
#' @param size Point size. Default = \code{3}.
#' @param return.DR.info Final data combining metadata and DR info should be returned. Default = \code{FALSE}.
#' @keywords dimensional reduction
#' @keywords median expression values
#' @keywords kmeans
#' @export median.dr
#' @import ggfortify
#' @importFrom SummarizedExperiment assay
#' @importFrom utils capture.output
#' @importFrom stats kmeans
#' @examples
#' \dontrun{
#' median.dr(fcs, color.by = "condition")
#' }

median.dr <- function(fcs.SCE, assay.i = "normalized", markers.to.use = "all", dr.method = "kmeans", num.k = 3, perplexity.tsne = 15, n.neighbors.umap = 50, color.by = "filename", label.by = NULL, size = 3, return.DR.info = F){
  data <- t(assay(fcs.SCE, i = assay.i))
  metadata <- fcs.SCE@metadata$reduced_metadata

  if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- colnames(data)

  ## prepare median tables
  med <- median.values(fcs.SCE, assay.i = assay.i)

  ## dr
  if(grepl("^pca$|^tsne$|^umap$", tolower(dr.method))){
    invisible(capture.output(dr <- dim.reduction(med, dr.method = dr.method, markers.to.use = markers.to.use,
                                                 perplexity.tsne = perplexity.tsne, n.neighbors.umap = n.neighbors.umap)))

    ## plotting
    drmd <- merge(metadata, dr[[1]], by = "row.names")[-1]
    dr.plotting(drmd, plot.dr = dr.method, color.by = color.by, label.by = label.by,
                size = size, raster = F)
  }else if(tolower(dr.method) == "kmeans"){
    dr <- kmeans(med, num.k)
    drmd <- merge(med, as.data.frame(dr$cluster), by = "row.names")[,-1]

    print(autoplot(dr, data = med, frame = T, frame.type = "norm") + theme_bw())
  }else{
    stop("Please, specify a valid reduction method: PCA, tSNE, UMAP or kmeans.", call. = F)
  }

  if(return.DR.info) return(drmd)
}
