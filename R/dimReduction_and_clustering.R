#' dim.reduction
#'
#' This function calculates cluster proportions for each identified cluster and draws it over a stacked barplot by condition.
#' @param expr_data A data.frame with marker expression values for each sample, patient, cell...
#' @param metadata A data.frame with additional information for \code{expr_data}.
#' @param reduction_method Dimensional reduction method. Possible values are "PCA", "tSNE", "UMAP" or "all".
#' @param num_threads Number of threads. Default = 4
#' @param tsne_perplexity Value for perplexity parameter in tSNE calculation (\href{https://distill.pub/2016/misread-tsne/}{more information}). Default = 100
#' @param set.seed \code{\link[base:set.seed]{set.seed()}} for dimensional reduction. Default = 333
#' @keywords dimensional reduction
#' @keywords tSNE
#' @keywords PCA
#' @keywords UMAP
#' @return The final output is a list with two elements \enumerate{
#'           \item \code{dr}: A data.frame with input table together with dimensional reduction data,
#'           \item \code{dr_melted}: A data.frame with input table together with dimensional reduction data but collapsed to long format by marker.}
#' @export dim.reduction
#' @method dim reduction
#' @examples
#' \dontrun{
#' dr <- dim.reduction(sel_expr, metadata = sel_exprL, reduction_method = "all",
#'     set.seed = 1234)
#'     }

dim.reduction <- function(expr_data, metadata, reduction_method = c("all", "PCA", "tSNE", "UMAP"), num_threads = 4, tsne_perplexity = 100, set.seed = 333){
  set.seed(set.seed)

  if(reduction_method == "PCA"){
    PCA_out <- stats::prcomp(expr_data, center = TRUE, scale. = FALSE)

    #dataframe for visualization
    dr <- data.frame(metadata, expr_data, PCA1 = PCA_out$x[,1], PCA2 = PCA_out$x[,2])
    drm <- data.table::melt(dr, measure.vars = colnames(expr_data), value.name = "expression", variable.name = "antigen")

  }else if(reduction_method == "tSNE"){
    tsne_out <- Rtsne.multicore::Rtsne.multicore(expr_data, num_threads = num_threads, check_duplicates = FALSE, pca = FALSE, perplexity = tsne_perplexity)

    #dataframe for visualization
    dr <- data.frame(metadata, expr_data, tSNE1 = tsne_out$Y[,1], tSNE2 = tsne_out$Y[,2])
    drm <- data.table::melt(dr, measure.vars = colnames(expr_data), value.name = "expression", variable.name = "antigen")

  }else if(reduction_method == "UMAP"){
    umap_out <- uwot::tumap(expr_data, n_neighbors = 50, init = "random", n_threads = num_threads, verbose = F)

    #dataframe for visualization
    dr <- data.frame(metadata, expr_data, UMAP1 = umap_out[,1], UMAP2 = umap_out[,2])
    drm <- data.table::melt(dr, measure.vars = colnames(expr_data), value.name = "expression", variable.name = "antigen")

  }else if(reduction_method == "all"){ 
    
    print("Calculating PCA reduction...")
    PCA_out <- stats::prcomp(expr_data, center = TRUE, scale. = FALSE)
    print("Calculating tSNE reduction...")
    tsne_out <- Rtsne.multicore::Rtsne.multicore(expr_data, num_threads = num_threads, check_duplicates = FALSE, pca = FALSE, perplexity = tsne_perplexity)
    print("Calculating UMAP reduction...")
    umap_out <- uwot::tumap(expr_data, n_neighbors = 50, init = "random", n_threads = num_threads, verbose = F)

    #dataframe for visualization
    dr <- data.frame(metadata, expr_data,
                     PCA1 = PCA_out$x[,1], PCA2 = PCA_out$x[,2],
                     tSNE1 = tsne_out$Y[,1], tSNE2 = tsne_out$Y[,2],
                     UMAP1 = umap_out[,1], UMAP2 = umap_out[,2])
    drm <- data.table::melt(dr, measure.vars = colnames(expr_data), value.name = "expression", variable.name = "antigen")

  }else{
    print('You must to provide a reduction method: PCA, tSNE, UMAP or all')
  }
  return(list(dr = dr, dr_melted = drm))
}


#' fsom.clustering
#'
#' This function calculates a SOM (Self-Organizing Map) clustering from a data.frame with marker expression values for all cells.
#' For additional information go to \href{https://github.com/SofieVG/FlowSOM}{\code{\link{FlowSOM}}} package.
#' @param data A data.frame with the marker expression values for each cell.
#' @param markers_to_use Markers to use in the clustering calculation. Default = "surface_markers"
#' @param markers_to_plot Markers to plot in the final image. Default = "tree" (it plots a simple tree, without markers).
#' @param set.seed \code{\link[base:set.seed]{set.seed()}} for clustering calculation. Default = 333
#' @param saveRDS Logical whether saving the final calculated object. Default = \code{FALSE}
#' @keywords SOM clustering
#' @export
#' @examples
#' \dontrun{
#' fsom <- fsom.clustering(fcs, markers_to_use = "surface_markers",
#'     markers_to_plot = "tree")
#' fsom <- fsom.clustering(cell_expressions, markers_to_use = "surface_markers",
#'     markers_to_plot = "CXCR3")
#'     }

fsom.clustering <- function(data, markers_to_use = "surface_markers",
                            markers_to_plot = "tree", set.seed = 333, saveRDS = F){
  set.seed(set.seed)

  if(class(data)[1] != "flowSet"){
    data <- premessa::as_flowFrame(as.matrix(data))
  }

  print("Calculate SOM...")
  fsom <- FlowSOM::ReadInput(data, transform = FALSE, scale = FALSE) #read data
  fsom <- FlowSOM::BuildSOM(fsom, colsToUse = eval(parse(text = markers_to_use))) #build SOM

  if(saveRDS){
    saveRDS(fsom, "som.rds", compress = FALSE)
  }

  print("Build MST...")
  fsom <- FlowSOM::BuildMST(fsom, tSNE = TRUE, silent = T) #build MST for visualization of clustering

  #plot the MST to evaluate the marker fluorescence (or general tree) intensity for each SOM
  if(markers_to_plot == "tree"){
    FlowSOM::PlotStars(fsom)
  }else{
    for(marker in markers_to_plot) FlowSOM::PlotMarker(fsom, marker)
  }

  return(fsom)
}


#' fsom.metaclustering
#'
#' This function calculates a consensus clustering to reduce the number of previous identified clusters through SOM (Self-Organizing Map) clustering.
#' For additional information go to \href{http://bioconductor.org/packages/release/bioc/vignettes/ConsensusClusterPlus/inst/doc/ConsensusClusterPlus.pdf}{\code{\link{ConsensusClusterPlus}}} package.
#' @param fsom Object previously calculated with \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}}.
#' @param num_clusters_metaclustering Number of clusters to reduce. Default = 40
#' @param plotting Logical whether plotting final SOM recalculation. Default = \code{TRUE}
#' @param folder_name Name for created folder where consensus clustering will be stored. Default = "consensus_plot"
#' @param set.seed \code{\link[base:set.seed]{set.seed()}} for clustering calculation. Default = 333
#' @param saveRDS Logical whether saving final calculated object. Default = \code{FALSE}
#' @keywords SOM clustering
#' @keywords metaclustering
#' @keywords consensus clustering
#' @export
#' @examples
#' \dontrun{metaclusters <- fsom.metaclustering(fsom = fsom, num_clusters_metaclustering = 40)}

fsom.metaclustering <- function(fsom, num_clusters_metaclustering = 40,
                                plotting = T, folder_name = "consensus_plots",
                                set.seed = 333, saveRDS = F){

  print("Metaclustering...")
  mc <- ConsensusClusterPlus::ConsensusClusterPlus(t(fsom$map$codes), maxK = num_clusters_metaclustering, reps = 100,
                             pItem = 0.9, pFeature = 1, title = folder_name, plot = "png",
                             clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = set.seed, verbose = F)

  #get cluster ids for each cell
  code_clustering1 <- mc[[num_clusters_metaclustering]]$consensusClass %>% as.factor()
  cell_clustering1 <- code_clustering1[fsom$map$mapping[,1]]

  #add clustering to original MST and color by cluster colors
  if(plotting){
    FlowSOM::PlotStars(fsom, backgroundValues = code_clustering1, backgroundColor = alpha(colors_palette, alpha = 0.4))
  }

  return(list(metaclusters = cell_clustering1, plotStars_value = code_clustering1))
}
