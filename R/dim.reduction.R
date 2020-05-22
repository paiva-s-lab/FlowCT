#' dim.reduction
#'
#' It calculates a dimensional reduction (DR) from a \code{fcs.SCE} object (or an expression table). Three different DR methods are available: Principal Component Analsis (PCA), t-Distributed Stochastic Neighbor Embedding (t-SNE) and  Uniform Manifold Approximation and Projection (UMAP).
#' @param data A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}} or a expression table with events in rows and markers in columns.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate the DR (this option is useless if input is not a \code{fcs.SCE} object). Default = \code{"normalized"}.
#' @param markers.to.use Markers to take account in the DR calculus. Default = \code{"all"}.
#' @param dr.method DR method to calculate. Possible values are "PCA", "tSNE" and/or "UMAP".
#' @param num.threads Number of threads for DR calculus. If \code{NULL} (default), all cores available minus one will be used.
#' @param perplexity.tsne Value for perplexity parameter in tSNE calculation (\href{https://distill.pub/2016/misread-tsne/}{more information}). Default = \code{100}.
#' @param n.neighbors.umap Value for neighbors parameter in UMAP calculation (\href{https://www.math.upenn.edu/~jhansen/2018/05/04/UMAP/}{more information}). Default = \code{50}.
#' @keywords dimensional reduction
#' @keywords tSNE
#' @keywords PCA
#' @keywords UMAP
#' @importFrom SummarizedExperiment colData assay
#' @importFrom parallel detectCores
#' @importFrom stats prcomp
#' @importFrom Rtsne Rtsne
#' @importFrom uwot tumap
#' @importFrom SingleCellExperiment reducedDims
#' @import dplyr
#' @export dim.reduction
#' @method dim reduction
#' @examples
#' \dontrun{
#' fcs <- dim.reduction(fcs, dr.method = "tSNE", markers.to.use = c("CD8", "CD27", "CCR4", "CD45RA", "CD4"))
#' fcs <- dim.reduction(fcs, dr.method = "UMAP", n.neighbors.umap = 10)
#' fcs <- dim.reduction(fcs, dr.method = c("pca", "Umap"))
#' }

dim.reduction <- function(data, assay.i = "normalized", markers.to.use = "all", dr.method, num.threads = NULL, perplexity.tsne = 100, n.neighbors.umap = 50){
  if(!all(grepl("^pca$|^tsne$|^umap$", tolower(dr.method)))) stop('The available reduction methods are: PCA, tSNE and/or UMAP.\n', call. = F)
  
  if(is.null(num.threads)) num.threads <- detectCores()-1
  
  if(class(data) == "SingleCellExperiment"){
    if(markers.to.use == "all") markers.to.use <- rownames(data)
    data1 <- t(assay(data, i = assay.i))[,markers.to.use]
  }else{
    if(markers.to.use == "all") markers.to.use <- colnames(data)
    data1 <- data[,markers.to.use]
  }
  
  drs <- list()
  if("pca" %in% tolower(dr.method)){
    cat(">>> PCA calculation...\n")
    drs[["PCA"]] <- prcomp(data1, center = TRUE, scale. = FALSE)$x[,1:2]
    colnames(drs[["PCA"]]) <- paste0("pca", c(1:2))
  }
  if("tsne" %in% tolower(dr.method)){
    cat(">>> tSNE calculation...\n")
    drs[["tSNE"]] <- Rtsne(data1, num_threads = num.threads, check_duplicates = FALSE, pca = FALSE, perplexity.tsne = perplexity.tsne)$Y
    colnames(drs[["tSNE"]]) <- paste0("tsne", c(1:2))
  } 
  if("umap" %in% tolower(dr.method)){
    cat(">>> UMAP calculation...\n")
    drs[["UMAP"]] <- tumap(data1, n_neighbors = n.neighbors.umap, init = "random", n_threads = num.threads, verbose = F)
    colnames(drs[["UMAP"]]) <- paste0("umap", c(1:2))
  }
  
  if(class(data) == "SingleCellExperiment"){
    reducedDims(data) <- drs
    return(data)
  }else{
    return(drs)
  }
}
