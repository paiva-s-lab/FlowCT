#' Calculate a DR with median values
#'
#' This function draws a plot from internally calculated K-means and dimensional reduction (DR).
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param dr.method DR method to calculate. Possible values are "PCA", "tSNE" or "UMAP". Take into account the number of samples for this DR, for few samples t-SNE and UMAP are useless. Default = \code{"PCA"}.
#' @param num.k Number of clusters if \code{dr.method = kmeans}. Default = \code{3}.
#' @param kmeans.frame Should kmeans cluster be sorrunded by a circular frame within the DR plot?. Default = \code{TRUE}.
#' @param perplexity.tsne Value for perplexity parameter in tSNE calculation (\href{https://distill.pub/2016/misread-tsne/}{more information}). For little number of samples this DR (t-SNE) cannot be calculated because preplexity restrictions. Default = \code{15}.
#' @param n.neighbors.umap Value for neighbors parameter in UMAP calculation (\href{https://www.math.upenn.edu/~jhansen/2018/05/04/UMAP/}{more information}). Default = \code{50}.
#' @param color.by Variable name (from \code{colData(fcs.SCE)}) for dots coloring. Default = \code{"filename"}.
#' @param shape.by Variable name (from \code{colData(fcs.SCE)}) for dots shaping. Default = \code{"NULL"}.
#' @param label.by Variable name (from \code{colData(fcs.SCE)}) for labelling. Default = \code{NULL}.
#' @param size Point size. Default = \code{3}.
#' @param return.DR.info Final data combining metadata, K-means and DR info should be returned?. Default = \code{FALSE}.
#' @keywords dimensional reduction
#' @keywords median expression values
#' @keywords kmeans
#' @export median.dr
#' @import ggplot2
#' @importFrom SummarizedExperiment assay
#' @importFrom utils capture.output
#' @importFrom stats kmeans
#' @examples
#' \dontrun{
#' median.dr(fcs, color.by = "condition")
#' }

median.dr <- function(fcs.SCE, assay.i = "normalized", markers.to.use = "all", dr.method = "PCA", num.k = 3, kmeans.frame = T, perplexity.tsne = 15, n.neighbors.umap = 50, color.by = "filename", shape.by = NULL, label.by = NULL, size = 3, return.DR.info = F){
  data <- t(assay(fcs.SCE, i = assay.i))
  metadata <- fcs.SCE@metadata$reduced_metadata

  if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- colnames(data)

  ## prepare median tables and k-means
  med <- median.values(fcs.SCE, assay.i = assay.i)
  set.seed(333); k <- kmeans(med, num.k)
  mdk <- merge(metadata, as.data.frame(as.factor(k$cluster)), by = "row.names")[,-1]
  colnames(mdk)[ncol(mdk)] <- paste0("kmeans", ".k", num.k)
  rownames(mdk) <- mdk$filename

  ## dr
  if(grepl("^pca$|^tsne$|^umap$", tolower(dr.method))){
    invisible(capture.output(dr <- dim.reduction(med, dr.method = dr.method, markers.to.use = markers.to.use,
                                                 perplexity.tsne = perplexity.tsne, n.neighbors.umap = n.neighbors.umap)))
    colnames(dr[[1]]) <- paste0("dr", 1:2)

    ## plotting
    # drmd <- merge(mdk, dr[[1]], by = "row.names")[-1]
    drmd <- cbind(mdk, dr[[1]])
    if(!kmeans.frame){
      g <- dr.plotting(drmd, plot.dr = dr.method, color.by = color.by, shape.by = shape.by, label.by = label.by,
                  size = size, raster = F)
    }else{
      g <- ggplot(drmd, aes_string("dr1", "dr2", color = paste0("kmeans", ".k", num.k), shape = shape.by)) +
              geom_point() +
              stat_ellipse(geom = "polygon", alpha = 1/4, type = "t", level = 0.8,
                           inherit.aes = F, aes_string("dr1", "dr2", fill = paste0("kmeans", ".k", num.k),
                                                       color = paste0("kmeans", ".k", num.k))) +
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA))
    }
  }else{
    stop("Please, specify a valid reduction method: PCA, tSNE, UMAP or kmeans.", call. = F)
  }

  if(return.DR.info){
    colnames(drmd)[grepl("dr", colnames(drmd))] <- paste0(dr.method, 1:length(grep("dr", colnames(drmd))))
    return(list(data = drmd, plot = g))
  }else{
    return(g)
  }
}
