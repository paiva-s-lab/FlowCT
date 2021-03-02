#' Calculate a DR with median values
#'
#' This function performs a dimensional reduction over median expressino values for each file within the \code{fcs.SCE} object and, if user specifies, calculated K-means clustering.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param markers.to.use Vector with markers to use. By default (\code{"all"}), all markers will be used.
#' @param dr.method DR method to calculate. Possible values are "PCA" (default), "tSNE", "UMAP", "DENSMAP" or "DENSNE". Take into account the number of samples for this DR, for few samples t-SNE and UMAP are useless. If `pca.loadings`, the plotting result will show the weight of each marker for each PC (the first two components).
#' @param num.k Number of clusters if \code{dr.method = kmeans}. If \code{NULL} (default), a classical plot is drawn, but if a number is indicated, K-means clustering is performed and displayed as ellipse-rounded groups.
#' @param perplexity.tsne Value for perplexity parameter in tSNE calculation (\href{https://distill.pub/2016/misread-tsne/}{more information}). For little number of samples this DR (t-SNE) cannot be calculated because preplexity restrictions. Default = \code{15}.
#' @param n.neighbors.umap Value for neighbors parameter in UMAP calculation (\href{https://www.math.upenn.edu/~jhansen/2018/05/04/UMAP/}{more information}). Default = \code{50}.
#' @param color.by Variable name (from \code{colData(fcs.SCE)}) for dots coloring. Default = \code{"filename"}.
#' @param shape.by Variable name (from \code{colData(fcs.SCE)}) for dots shaping. Default = \code{"NULL"}.
#' @param label.by Variable name (from \code{colData(fcs.SCE)}) for labelling. Default = \code{NULL}.
#' @param size Point size. Default = \code{3}.
#' @param return.DR Final data combining metadata, K-means and DR info should be returned?. Default = \code{FALSE}.
#' @keywords dimensional reduction
#' @keywords median expression values
#' @keywords kmeans
#' @export median.dr
#' @import ggplot2
#' @importFrom SummarizedExperiment assay
#' @importFrom stats kmeans
#' @importFrom utils capture.output
#' @examples
#' \dontrun{
#' median.dr(fcs, color.by = "condition")
#' }

median.dr <- function(fcs.SCE, assay.i = "normalized", markers.to.use = "all", dr.method = "PCA", num.k = NULL, perplexity.tsne = 15, n.neighbors.umap = 50, color.by = "filename", shape.by = NULL, label.by = NULL, size = 3, return.DR = F){
  if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- rownames(fcs.SCE)
  metadata <- dplyr::distinct(as.data.frame(colData(fcs.SCE)), .data$filename, .keep_all = T)
  rownames(metadata) <- metadata$filename
  
  med <- median.values(fcs.SCE, assay.i = assay.i)

  if(grepl("^pca$|^tsne$|^umap$|^densmap$|^densne$", tolower(dr.method))){
    invisible(capture.output(dr <- dim.reduction(med, dr.method = dr.method, markers.to.use = markers.to.use,
                                                 perplexity.tsne = perplexity.tsne, n.neighbors.umap = n.neighbors.umap)))
    drmd <- cbind(metadata, dr$DRs[[1]])
    
    ## plotting
    if(is.null(num.k)){
      g <- dr.plotting(drmd, plot.dr = dr.method, color.by = color.by, shape.by = shape.by, label.by = label.by,
                       size = size)
    }else{
    set.seed(333); k <- kmeans(med, num.k)
    drmd <- merge(drmd, as.data.frame(as.factor(k$cluster)), by.x = "filename", by.y = "row.names")
    colnames(drmd)[ncol(drmd)] <- paste0("kmeans", ".k", num.k)

    g <- dr.plotting(drmd, plot.dr = dr.method, color.by = color.by, shape.by = shape.by, label.by = label.by, size = size) + 
              stat_ellipse(geom = "polygon", alpha = 1/4, type = "t", level = 0.8,
                           inherit.aes = F, aes_string(colnames(dr[[1]])[1], colnames(dr[[1]])[2], 
                                          fill = paste0("kmeans", ".k", num.k)))


    }
  }else if(tolower(dr.method) == "pca.loadings"){
    invisible(capture.output(dr <- dim.reduction(med, dr.method = "pca", markers.to.use = markers.to.use)))
    return(dr.plotting(dr, plot.dr = "pca.loadings"))
  }else{
    stop("Please, specify a valid reduction method: PCA, tSNE, UMAP, DENSMAP or DENSNE.", call. = F)
  }
  
  if(return.DR){
    colnames(drmd)[grepl("dr", colnames(drmd))] <- paste0(dr.method, 1:length(grep("dr", colnames(drmd))))
    return(list(data = drmd, plot = g))
  }else{
    return(g)
  }
}
