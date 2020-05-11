#' dim.reduction
#'
#' It calculates a dimensional reduction (DR) from a \code{FCS.SE} object (or an expression table). Three different DR methods are available: Principal Component Analsis (PCA), t-Distributed Stochastic Neighbor Embedding (t-SNE) and  Uniform Manifold Approximation and Projection (UMAP).
#' @param data A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}} or a expression table with events in rows and markers in columns.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate the DR (this option is useless if input is not a \code{FCS.SE} object). Default = \code{"normalized"}.
#' @param metadata If \code{data} is an expression matrix, an additional metadata table with information for each FCS file is necessary.
#' @param markers.to.use Markers to take account in the DR calculus. Default = \code{"all"}.
#' @param dr.method DR method to calculate. Possible values are "PCA", "tSNE" or "UMAP", individually, or "all".
#' @param num.threads Number of threads for DR calculus. If \code{NULL} (default), all cores available minus one will be used.
#' @param perplexity.tsne Value for perplexity parameter in tSNE calculation (\href{https://distill.pub/2016/misread-tsne/}{more information}). Default = \code{100}.
#' @param n.neighbors.umap Value for neighbors parameter in UMAP calculation (\href{https://www.math.upenn.edu/~jhansen/2018/05/04/UMAP/}{more information}). Default = \code{50}.
#' @keywords dimensional reduction
#' @keywords tSNE
#' @keywords PCA
#' @keywords UMAP
#' @return The final output is an object of type \code{list} with two elements \enumerate{
#'           \item \code{dr}: A data.frame combining metadata and expression values, together with dimensional reduction data.
#'           \item \code{dr_melted}: A data.frame combining metadata and expression values, together with dimensional reduction data but transformed to long format by marker for later plotting.}
#' @importFrom SummarizedExperiment colData assay
#' @export dim.reduction
#' @method dim reduction
#' @examples
#' \dontrun{
#' dr_tsne <- dim.reduction(fcs_se, dr.method = "tSNE", markers.to.use = c("CD8", "CD27", "CCR4", "CD45RA", "CD4"))
#' dr_umap <- dim.reduction(fcs_se, dr.method = "UMAP", n.neighbors.umap = 10)
#' dr_all <- dim.reduction(fcs_se, dr.method = "all")
#' }

dim.reduction <- function(data, assay.i = "normalized", metadata = NULL, markers.to.use = "all", dr.method = c("all", "PCA", "tSNE", "UMAP"), num.threads = NULL, perplexity.tsne = 100, n.neighbors.umap = 50){
  if(is.null(num.threads)) num.threads <- parallel::detectCores()-1
  if(markers.to.use == "all") markers.to.use <- rownames(data)

  if(class(data) == "SummarizedExperiment"){
    metadata <- colData(data)
    data <- t(assay(data, i = assay.i))
  }else{
    metadata <- metadata
  }
  
  ## DR
  if(dr.method == "PCA"){
    PCA_out <- stats::prcomp(data[,markers.to.use], center = TRUE, scale. = FALSE)
    
    #dataframe for visualization
    dr <- data.frame(metadata, data[,markers.to.use], PCA1 = PCA_out$x[,1], PCA2 = PCA_out$x[,2])
    drm <- data.table::melt(dr, measure.vars = colnames(data[,markers.to.use]), value.name = "expression", variable.name = "antigen")
    
  }else if(dr.method == "tSNE"){
    tsne_out <- Rtsne::Rtsne(data[,markers.to.use], num_threads = num.threads, check_duplicates = FALSE, pca = FALSE, perplexity.tsne = perplexity.tsne)
    
    #dataframe for visualization
    dr <- data.frame(metadata, data[,markers.to.use], tSNE1 = tsne_out$Y[,1], tSNE2 = tsne_out$Y[,2])
    drm <- data.table::melt(dr, measure.vars = colnames(data[,markers.to.use]), value.name = "expression", variable.name = "antigen")
    
  }else if(dr.method == "UMAP"){
    umap_out <- uwot::tumap(data[,markers.to.use], n_neighbors = n.neighbors.umap, init = "random", n_threads = num.threads, verbose = F)
    
    #dataframe for visualization
    dr <- data.frame(metadata, data[,markers.to.use], UMAP1 = umap_out[,1], UMAP2 = umap_out[,2])
    drm <- medata.table::melt(dr, measure.vars = colnames(data[,markers.to.use]), value.name = "expression", variable.name = "antigen")
    
  }else if(dr.method == "all"){
    cat("Calculating PCA reduction...\n\n")
    PCA_out <- stats::prcomp(data[,markers.to.use], center = TRUE, scale. = FALSE)
    
    cat("Calculating tSNE reduction...\n\n")
    tsne_out <- Rtsne::Rtsne(data[,markers.to.use], num_threads = num.threads, check_duplicates = FALSE, pca = FALSE, perplexity.tsne = perplexity.tsne)
    
    cat("Calculating UMAP reduction...\n\n")
    umap_out <- uwot::tumap(data[,markers.to.use], n_neighbors = n.neighbors.umap, init = "random", n_threads = num.threads, verbose = F)
    
    #dataframe for visualization
    dr <- data.frame(metadata, data[,markers.to.use], 
                     PCA1 = PCA_out$x[,1], PCA2 = PCA_out$x[,2],
                     tSNE1 = tsne_out$Y[,1], tSNE2 = tsne_out$Y[,2],
                     UMAP1 = umap_out[,1], UMAP2 = umap_out[,2])
    drm <- data.table::melt(dr, measure.vars = colnames(data[,markers.to.use]), value.name = "expression", variable.name = "antigen")
    
  }else{
    cat('>>> You must to provide a reduction method: PCA, tSNE, UMAP or all\n')
  }
  
  return(list(dr = dr, dr_melted = drm))
}
