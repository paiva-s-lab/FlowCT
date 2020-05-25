#' fsom.clustering
#'
#' This function calculates a Self-Organizing Map (SOM) clustering and a simultaneous metaclustering from a \code{fcs.SCE} object. For additional information go to \href{https://github.com/SofieVG/FlowSOM}{\code{FlowSOM} package}.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate SOM clustering. Default = \code{"normalized"}.
#' @param markers.to.use Markers used for the clustering calculation. Default = \code{"all"}.
#' @param markers.to.plot Markers to plot in the final minimum spanning tree (MST). Posible values are \code{NULL} (no plotting), \code{"tree"} (plotting a general MST with all identified clusters), \code{"tree_metaclustering"} (plotting a general MST but colored by metaclustering results) or a vector with markers (and it draws multiple MSTs for each marker indicated).
#' @param k.metaclustering Number of metaclusters to calculate in the metaclustering step. Default = \code{40}.
#' @param metaclustering.name Name for PDF with metaclustering results (see \href{https://bioconductor.org/packages/release/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf}{\code{ConsensusClusterPlus} package} for interpretation help). Default = \code{NULL}.
#' @keywords SOM clustering
#' @keywords SOM metaclustering
#' @keywords minimum spanning tree MST
#' @export
#' @import FlowSOM
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @examples
#' \dontrun{
#' fsom <- fsom.clustering(fcs.SCE = fcs, markers.to.use = c("CD8", "CD27", "CCR4", "CD45RA", "CD4"), k.metaclustering = 20)
#' fcs$SOM <- fsom$metaclusters #add SOM information to colData(fcs) as a new column
#' }

fsom.clustering <- function(fcs.SCE, assay.i = "normalized", markers.to.use = "all", markers.to.plot = NULL, k.metaclustering = NULL, metaclustering.name = NULL){
  set.seed(333)
  
  data <- as.flowSet.SE(fcs.SCE, assay.i)
  if(markers.to.use == "all") markers.to.use2 <- rownames(fcs.SCE) else markers.to.use2 <- markers.to.use
  
  ## FSOM clustering
  cat("Calculating SOM clustering...\n")
  fsom <- suppressMessages(ReadInput(data, transform = FALSE, scale = FALSE)) #read data
  fsom <- suppressMessages(BuildSOM(fsom, colsToUse = markers.to.use2)) #build SOM
  
  cat("Building MST...\n")
  fsom <- suppressMessages(BuildMST(fsom, tSNE = TRUE, silent = T)) #build MST for visualization of clustering
  
  ## Metaclustering
  if(!is.null(k.metaclustering)){
    cat("Metaclustering...\n")
    if(!is.null(metaclustering.name)){
      mc <- suppressMessages(ConsensusClusterPlus(t(fsom$map$codes), maxK = k.metaclustering, reps = 100,
                                                  pItem = 0.9, pFeature = 1, title = metaclustering.name, plot = "pdf",
                                                  clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                                                  distance = "euclidean", seed = 333, verbose = F))
    }else{
      mc <- suppressMessages(ConsensusClusterPlus(t(fsom$map$codes), maxK = k.metaclustering, reps = 100,
                                                  pItem = 0.9, pFeature = 1, title = "consensus_plots", plot = "pdf",
                                                  clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                                                  distance = "euclidean", seed = 333, verbose = F))
      unlink("consensus_plots", recursive = TRUE)
    }
    
    #get cluster ids for each cell
    code_clustering1 <- mc[[k.metaclustering]]$consensusClass %>% as.factor()
    cell_clustering1 <- code_clustering1[fsom$map$mapping[,1]]
    
    ## plotting
    #plot the MST to evaluate the marker fluorescence (or general tree) intensity for each SOM
    if(!is.null(markers.to.plot)){
      if(markers.to.plot == "tree"){
        PlotStars(fsom)
      }else if(markers.to.plot == "tree_metaclustering"){
        PlotStars(fsom, backgroundValues = code_clustering1, backgroundColor = alpha(div.colors(length(code_clustering1)), alpha = 0.7))
      }else{
        for(marker in markers.to.plot) PlotMarker(fsom, marker)
      }
    }
    
    return(list(fsom = fsom, metaclusters = cell_clustering1, plotStars_value = code_clustering1))
  }else{
    return(fsom)
  }
}
