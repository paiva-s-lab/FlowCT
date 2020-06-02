#' Calculate dimensional reductions
#'
#' It calculates a dimensional reduction (DR) from a \code{fcs.SCE} object (or an expression table). Three different DR methods are available: Principal Component Analsis (PCA), t-Distributed Stochastic Neighbor Embedding (t-SNE) and  Uniform Manifold Approximation and Projection (UMAP).
#' @param data A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}} or a expression table with events in rows and markers in columns.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate the DR (this option is useless if input is not a \code{fcs.SCE} object). Default = \code{"normalized"}.
#' @param markers.to.use Markers to take account in the DR calculus. Default = \code{"all"}.
#' @param dr.method DR method to calculate. Possible values are "PCA", "tSNE" and/or "UMAP".
#' @param dr.suffix Suffix to add to DR calculation and avoid results overwritting. If \code{dr.suffix = "default"} (default option), an additional tag with DR internal parameters will be added to DR name (e.g., "tSNE.px100" would refer to tSNE calculation with a perplexity value of 100); if \code{dr.suffix = "none"}, only DR name will be kept, and if user specifies a suffix, it will be added to DR name in a similar way to "default" option.
#' @param num.threads Number of threads for DR calculus. If \code{NULL} (default), all cores available minus one will be used. If this option is enable, reproducibility could be comprised.
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

dim.reduction <- function(data, assay.i = "normalized", markers.to.use = "all", dr.method, dr.suffix = "default", num.threads = NULL, perplexity.tsne = 100, n.neighbors.umap = 50){
  if(!all(grepl("^pca$|^tsne$|^umap$", tolower(dr.method)))) stop('The available reduction methods are: PCA, tSNE and/or UMAP.\n', call. = F)
  if(is.null(num.threads)) num.threads <- detectCores()-1

  if(class(data)[1] == "SingleCellExperiment"){
    if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- rownames(data)
    data1 <- t(assay(data, i = assay.i))[,markers.to.use]
  }else{
    if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- colnames(data)
    data1 <- data[,markers.to.use]
  }

  drs <- list()
  if("pca" %in% tolower(dr.method)){
    cat(">>> PCA calculation...\n")
    set.seed(333); drs[["PCA"]] <- prcomp(data1, center = TRUE, scale. = FALSE)$x[,1:2]
    colnames(drs[["PCA"]]) <- paste0("pca", c(1:2))
  }
  if("tsne" %in% tolower(dr.method)){
    cat(">>> tSNE calculation...\n")
    if(dr.suffix == "default"){
      drname <- paste0("tSNE.px", perplexity.tsne)
    }else if(dr.suffix == "none"){
      drname <- "tSNE"
    }else{
      drname <- paste0("tSNE.", dr.suffix)
    }
    set.seed(333); drs[[drname]] <- Rtsne(data1, num_threads = num.threads, check_duplicates = FALSE, pca = FALSE, perplexity.tsne = perplexity.tsne)$Y
    colnames(drs[[drname]]) <- paste0("tsne", c(1:2))

  }
  if("umap" %in% tolower(dr.method)){
    cat(">>> UMAP calculation...\n")
    if(dr.suffix == "default"){
      drname <- paste0("UMAP.neig", n.neighbors.umap)
    }else if(dr.suffix == "none"){
      drname <- "UMAP"
    }else{
      drname <- paste0("UMAP.", dr.suffix)
    }
    set.seed(333); drs[[drname]] <- tumap(data1, n_neighbors = n.neighbors.umap, init = "random", n_threads = num.threads, verbose = F)
    colnames(drs[[drname]]) <- paste0("umap", c(1:2))
  }

  if(class(data)[1] == "SingleCellExperiment"){
    ## combine DRs if another DR has been calculated
    if(length(names(data@int_colData@listData$reducedDims)) != 0){
      if(data@int_colData@listData$reducedDims@metadata$DR_assay != assay.i){
        warning("You are using a different assay.i than previous calculated DR, results have been overwritten with this new assay.", call. = F, immediate. = T)
        reducedDims(data) <- NULL
        reducedDims(data) <- drs
        data@int_colData@listData$reducedDims@metadata$DR_assay <- assay.i #assay used for DR

        return(data)
      }
      isecDR <- intersect(names(drs), names(data@int_colData@listData$reducedDims))
      if(length(isecDR) == 0){
        drs <- c(as.list(data@int_colData@listData$reducedDims), drs)
      }else{
        drs <- c(as.list(data@int_colData@listData$reducedDims), drs[!grepl(paste(isecDR, collapse = "|"), names(drs))])
      }
    }
    reducedDims(data) <- drs
    data@int_colData@listData$reducedDims@metadata$DR_assay <- assay.i #assay used for DR

    return(data)
  }else{
    return(drs)
  }
}
