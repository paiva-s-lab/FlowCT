#' Calculate dimensional reductions
#'
#' It calculates a dimensional reduction (DR) from a \code{fcs.SCE} object (or an expression table). Five different DR methods are available: Principal Component Analsis (PCA), t-Distributed Stochastic Neighbor Embedding (t-SNE),Density-preserving t-SNE (DENSNE), Uniform Manifold Approximation and Projection (UMAP) and Density-preserving UMAP (DENSMAP).
#' @param data A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}} or a expression table with events in rows and markers in columns.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate the DR (this option is useless if input is not a \code{fcs.SCE} object). Default = \code{"normalized"}.
#' @param markers.to.use Markers to take account in the DR calculus. Default = \code{"all"}.
#' @param dr.method DR method to calculate. Possible values are "PCA", "tSNE", "DENSNE", "UMAP" and/or "DENSMAP".
#' @param num.threads Number of threads for DR calculus. If \code{NULL} (default), all cores available minus one will be used. If this option is enable, reproducibility could be comprised.
#' @param perplexity.tsne Value for perplexity parameter in tSNE calculation (\href{https://distill.pub/2016/misread-tsne/}{more information}). Default = \code{100}.
#' @param n.neighbors.umap Value for neighbors parameter in UMAP calculation (\href{https://www.math.upenn.edu/~jhansen/2018/05/04/UMAP/}{more information}). Default = \code{50}.
#' @keywords dimensional reduction
#' @keywords tSNE
#' @keywords PCA
#' @keywords UMAP
#' @keywords DENSMAP
#' @keywords DENSNE
#' @importFrom SingleCellExperiment colData reducedDims
#' @importFrom SummarizedExperiment assay
#' @importFrom stats prcomp
#' @importFrom parallel detectCores
#' @import dplyr
#' @export dim.reduction
#' @method dim reduction
#' @examples
#' \dontrun{
#' fcs <- dim.reduction(fcs, dr.method = "tSNE", 
#'    markers.to.use = c("CD8", "CD27", "CCR4", "CD45RA", "CD4"))
#' fcs <- dim.reduction(fcs, dr.method = "UMAP", n.neighbors.umap = 10)
#' fcs <- dim.reduction(fcs, dr.method = c("pca", "Umap"))
#' }

## SingleCellExperiment method...
setGeneric("reducedDims<-", getGeneric("reducedDims<-", package = "SingleCellExperiment"))

dim.reduction <- function(data, assay.i = "normalized", markers.to.use = "all", dr.method, num.threads = NULL, perplexity.tsne = 100, n.neighbors.umap = 50,dens_frac = 0.5, dens_lambda = 0.5){
  if(!all(grepl("^pca$|^tsne$|^umap$|^densmap$|^densne$", tolower(dr.method)))) stop('The available reduction methods are: PCA, tSNE, UMAP, denSNE and/or densMAP.\n', call. = F)

  if(is.null(num.threads)) num.threads <- detectCores()-1

  if(class(data)[1] == "SingleCellExperiment"){
    if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- rownames(data)
    data1 <- t(assay(data, i = assay.i))[,markers.to.use]
  }else{
    if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- colnames(data)
    data1 <- data[,markers.to.use]
  }

  #drs <- list()
  drs<-reducedDims(data)
  if("pca" %in% tolower(dr.method)){
    cat(">>> PCA calculation...\n")
    drs[["PCA"]] <- prcomp(data1, center = TRUE, scale. = FALSE)$x[,1:2]
    colnames(drs[["PCA"]]) <- paste0("pca", c(1:2))
  }
  if("tsne" %in% tolower(dr.method)){
    if (!requireNamespace("Rtsne", quietly = TRUE)) {
      stop("Package \"Rtsne\" needed for this function to work. Please install it.", call. = FALSE)
    }

    cat(">>> tSNE calculation...\n")
    set.seed(333); drs[["tSNE"]] <- Rtsne::Rtsne(data1, num_threads = num.threads, check_duplicates = FALSE, pca = FALSE, perplexity.tsne = perplexity.tsne)$Y
    colnames(drs[["tSNE"]]) <- paste0("tsne", c(1:2))
  }
  if("umap" %in% tolower(dr.method)){
    if (!requireNamespace("uwot", quietly = TRUE)) {
      stop("Package \"uwot\" needed for this function to work. Please install it.", call. = FALSE)
    }

    cat(">>> UMAP calculation...\n")
    set.seed(333); drs[["UMAP"]] <- uwot::tumap(data1, n_neighbors = n.neighbors.umap, init = "random", n_threads = num.threads, verbose = F)
    colnames(drs[["UMAP"]]) <- paste0("umap", c(1:2))
  }
 if("densmap" %in% tolower(dr.method)){
    if (!requireNamespace("densvis", quietly = TRUE)) {
      stop("Package \"densvis\" needed for this function to work. Please install it.", call. = FALSE)
    }
    
    cat(">>> densmap calculation...\n")
    set.seed(333); drs[["DENSMAP"]] <- densvis::densmap(data1, dens_frac = dens_frac, dens_lambda = dens_lambda)
    colnames(drs[["DENSMAP"]]) <- paste0("densmap", c(1:2))
  }
  if("densne" %in% tolower(dr.method)){
    if (!requireNamespace("densvis", quietly = TRUE)) {
      stop("Package \"densvis\" needed for this function to work. Please install it.", call. = FALSE)
    }
    
    cat(">>> densne calculation...\n")
    set.seed(333); drs[["DENSNE"]] <- densvis::densne(data1, dens_frac = dens_frac, dens_lambda = dens_lambda)
    colnames(drs[["DENSNE"]]) <- paste0("densne", c(1:2))
  }
  if(class(data)[1] == "SingleCellExperiment"){
    reducedDims(data) <- drs
    return(data)
  }else{
    return(drs)
  }
}
