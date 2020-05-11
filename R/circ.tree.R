#' circ.tree
#'
#' This function plots a circular dendrogram and a heatmap according a \code{FCS.SE}. Every leaf has a different colored point size regarding the cell type and the frequency for each cluster.
#' For additional information go to \href{https://guangchuangyu.github.io/software/ggtree/documentation/}{\code{\link{ggtree}}} package.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which drawing the circular tree. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param dist.method Distance method measurement to be used. Possible values are "euclidean" (default), "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param hclust.method Hierarchical clustering method to be used. Possible values are "average" (default), "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid".
#' @param nodes If \code{"display"} (default), nodes will be numered. If contains a numeric vector with node numbers, areas defined from these nodes will be differently colored.
#' @param open.angle Angle aperture circular layout. Default = \code{100}.
#' @param dendro.labels Logical whether adding labels to dendrogram. Default = \code{FALSE}.
#' @param scale.size Numerical value indicating how much scale points in the dendogram terminal nodes. Default = \code{10}.
#' @keywords circular tree
#' @keywords dendrogram
#' @keywords nodes
#' @keywords hierachical clustering
#' @export
#' @import dplyr
#' @import phytools
#' @import ggtree
#' @import ggplot2
#' @importFrom treeio isTip
#' @examples
#' \dontrun{
#' # step 1: display all node numbers to select how to coloring areas
#' circ.tree(fcs.SE = fcs_seL, cell.clusters = fcs_seL$SOM_L_named, nodes = "display")
#' 
#' # step 2: color areas indicating node numbers
#' circ.tree(fcs.SE = fcs_seL, cell.clusters = fcs_seL$SOM_L_named, nodes = c(17, 25))
#' }

circ.tree <- function(fcs.SE, assay.i = "normalized", cell.clusters, dist.method = "euclidean", hclust.method = "average", 
                      nodes = "display", open.angle = 100, dendro.labels = FALSE, scale.size = 10){
  exprs <- t(SummarizedExperiment::assay(fcs.SE, i = assay.i))
  exprs_01 <- exprs.saturate(exprs)
  colors_palette <- div.colors(length(cell.clusters))
  
  ## Circular hyerarchical clustering tree
  #median expression of each marker for each cell population
  expr_medianL <- data.frame(exprs, cell_clustering = cell.clusters) %>%
    group_by(.data$cell_clustering) %>% summarize_all(list(stats::median)) %>% as.data.frame(.data)
  expr01_medianL <- data.frame(exprs_01, cell_clustering = cell.clusters) %>%
    group_by(.data$cell_clustering) %>% summarize_all(list(stats::median)) %>% as.data.frame(.data)
  
  #calculate cluster frequencies (and annotation for heatmap)
  clustering_propL <- data.frame(node = 1:length(levels(cell.clusters)), prop.table(table(cell.clusters))*100)
  
  #hierarchical clustering on clusters
  dL <- stats::dist(expr_medianL, method = dist.method)
  cluster_rowsL <- stats::hclust(dL, method = hclust.method)
  expr_heatL <- as.matrix(expr01_medianL)
  rownames(expr_heatL) <- expr01_medianL$cell_clustering
  
  #hyerarchical tree building
  hca <- ape::as.phylo(cluster_rowsL)
  hca$tip.label <- rownames(expr_heatL)
  
  if(nodes == "display"){
    print(ggtree(hca, layout = "fan", branch.length = 1) + geom_text2(aes_string(subset =! "isTip", label = "node"), hjust = -.3) + geom_tiplab())
  }else{
    p1 <- ggtree(hca, layout = "fan", open.angle = open.angle, branch.length = 1) 
    
    p1 <- p1 %<+% clustering_propL #add dataframe for geom_point level
    
    for(i in nodes){
      p1 <- p1 + geom_hilight(node = i, fill = sample(colors_palette, 1), alpha = .6) +
        geom_point(aes_string(shape = "1", color = "cell_clusters", size = "Freq")) +
        scale_size_area(max_size = scale.size) +
        scale_color_manual(values = colors_palette) + 
        theme(legend.position = "right")
    }
    
    if(dendro.labels){
      p1 + geom_tiplab2(offset=0.1, align = F, size=3)
    }
    
    expr_tree_plot <- expr01_medianL[,-1] #add saturated expression to tree heatmap
    rownames(expr_tree_plot) <- rownames(expr_heatL)
    
    print(gheatmap(p1, expr_tree_plot, offset = 0.5, width = 1, font.size = 2, colnames_angle = 0, hjust = 0,
                   colnames_position = "top", high = "#b30000", low = "#fff7f3") +
            theme(legend.position = NULL))
  }
}
