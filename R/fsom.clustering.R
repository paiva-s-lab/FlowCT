#' fsom.clustering
#'
#' This function calculates a Self-Organizing Map (SOM) clustering and a simultaneous metaclustering from a \code{FCS.SE} object. For additional information go to \href{https://github.com/SofieVG/FlowSOM}{\code{\link{FlowSOM}}} package.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate SOM clustering. Default = \code{"normalized"}.
#' @param markers.to.use Markers used for the clustering calculation. Default = \code{"all"}.
#' @param markers.to.plot Markers to plot in the final minimum spanning tree (MST). Posible values are \code{NULL} (no plotting), \code{"tree"} (plotting a general MST with all identified clusters) or a vector with markers (and it draws multiple MSTs for each marker indicated).
#' @param k.metaclustering Number of metaclusters to calculate in the metaclustering step. Default = \code{40}.
#' @param metaclustering.name Name for PDF with metaclustering results (see (\ref{https://bioconductor.org/packages/release/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf}{ConsensusClusterPlus package}) for interpretation help). Default = \code{NULL}.
#' @keywords SOM clustering
#' @keywords SOM metaclustering
#' @keywords minimum spanning tree MST
#' @export
#' @import FlowSOM
#' @examples
#' \dontrun{
#' fsom <- fsom.clustering(fcs.SE = fcs_se, markers.to.use = c("CD8", "CD27", "CCR4", "CD45RA", "CD4"), k.metaclustering = 20)
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
