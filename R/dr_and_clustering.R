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
#' @export dim.reduction
#' @method dim reduction
#' @examples
#' \dontrun{
#' dr_tsne <- dim.reduction(fcs_se, dr.method = "tSNE", 
#'     markers.to.use = c("CD8", "CD27", "CCR4", "CD45RA", "CD4"))
#' dr_umap <- dim.reduction(fcs_se, dr.method = "UMAP", n.neighbors.umap = 10)
#' dr_all <- dim.reduction(fcs_se, dr.method = "all")
#' }

dim.reduction <- function(data, assay.i = "normalized", metadata = NULL, markers.to.use = "all", dr.method = c("all", "PCA", "tSNE", "UMAP"), num.threads = NULL, perplexity.tsne = 100, n.neighbors.umap = 50){
  if(is.null(num.threads)) num.threads <- parallel::detectCores()-1
  if(markers.to.use == "all") markers.to.use <- rownames(data)

  if(class(data) == "SummarizedExperiment"){
    metadata <- SummarizedExperiment::colData(data)
    data <- t(SummarizedExperiment::assay(data, i = assay.i))
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


#' fsom.clustering
#'
#' This function calculates a Self-Organizing Map (SOM) clustering and a simultaneous metaclustering from a \code{FCS.SE} object. For additional information go to \href{https://github.com/SofieVG/FlowSOM}{\code{\link{FlowSOM}}} package.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate SOM clustering. Default = \code{"normalized"}.
#' @param markers.to.use Markers used for the clustering calculation. Default = \code{"all"}.
#' @param markers.to.plot Markers to plot in the final minimum spanning tree (MST). Posible values are \code{NULL} (no plotting), \code{"tree"} (plotting a general MST with all identified clusters) or a vector with markers (and it draws multiple MSTs for each marker indicated).
#' @param k.metaclustering Number of metaclusters to calculate in the metaclustering step. Default = \code{40}.
#' @param metaclustering.name Name for PDF with metaclustering results (see (\href{https://bioconductor.org/packages/release/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf}{ConsensusClusterPlus package}) for interpretation help). Default = \code{NULL}.
#' @keywords SOM clustering
#' @keywords SOM metaclustering
#' @keywords minimum spanning tree MST
#' @export
#' @import FlowSOM
#' @examples
#' \dontrun{
#' fsom <- fsom.clustering(fcs.SE = fcs_se, 
#'     markers.to.use = c("CD8", "CD27", "CCR4", "CD45RA", "CD4"), k.metaclustering = 20)
#' fcs_se$SOM <- fsom$metaclusters #add SOM information to colData(fcs.se) as a new column
#' }

fsom.clustering <- function(fcs.SE, assay.i = "normalized", markers.to.use = "all", markers.to.plot = NULL, k.metaclustering = 40, metaclustering.name = NULL){
  require(FlowSOM)

  data <- as.flowSet.SE(fcs.SE, assay.i)
  if(markers.to.use == "all") markers.to.use <- rownames(fcs.SE)
  
  ## FSOM clustering
  cat("Calculating SOM clustering...\n")
  fsom <- suppressMessages(ReadInput(data, transform = FALSE, scale = FALSE)) #read data
  fsom <- suppressMessages(BuildSOM(fsom, colsToUse = markers.to.use)) #build SOM
  
  cat("Building MST...\n")
  fsom <- suppressMessages(BuildMST(fsom, tSNE = TRUE, silent = T)) #build MST for visualization of clustering
  
  #plot the MST to evaluate the marker fluorescence (or general tree) intensity for each SOM
  if(!is.null(markers.to.plot)){
    if(markers.to.plot == "tree"){
      PlotStars(fsom)
    }else{
      for(marker in markers.to.plot) PlotMarker(fsom, marker)
    }
  }
  
  ## Metaclustering
  if(!is.null(k.metaclustering)){
    cat("Metaclustering...\n")
    if(!is.null(metaclustering.name)){
      mc <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(t(fsom$map$codes), maxK = k.metaclustering, reps = 100,
                                                  pItem = 0.9, pFeature = 1, title = metaclustering.name, plot = "pdf",
                                                  clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                                                  distance = "euclidean", seed = 333, verbose = F))
    }else{
      mc <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(t(fsom$map$codes), maxK = k.metaclustering, reps = 100,
                                                  pItem = 0.9, pFeature = 1, title = "consensus_plots", plot = "pdf",
                                                  clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                                                  distance = "euclidean", seed = 333, verbose = F))
      unlink("consensus_plots", recursive = TRUE)
    }
    
    #get cluster ids for each cell
    code_clustering1 <- mc[[k.metaclustering]]$consensusClass %>% as.factor()
    cell_clustering1 <- code_clustering1[fsom$map$mapping[,1]]
    
    #add clustering to original MST and color by cluster colors
    PlotStars(fsom, backgroundValues = code_clustering1, backgroundColor = alpha(div.colors(length(code_clustering1)), alpha = 0.7))
    
    return(list(fsom = fsom, metaclusters = cell_clustering1, plotStars_value = code_clustering1))
  }else{
    return(fsom)
  }
}
