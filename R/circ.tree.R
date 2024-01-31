#' Circular dendrogram tree
#'
#' This function plots a circular dendrogram and a heatmap according a \code{fcs.SCE}. Every leaf has a different colored point size regarding the cell type and the frequency for each cluster.
#' For additional information go to \href{https://guangchuangyu.github.io/software/ggtree/documentation/}{\code{ggtree} package}.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which drawing the circular tree. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param dist.method Distance method measurement to be used. Possible values are "euclidean" (default), "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param hclust.method Hierarchical clustering method to be used. Possible values are "average" (default), "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid".
#' @param nodes If \code{"display"} (default), nodes will be numered. If contains a numeric vector with node numbers, areas defined from these nodes will be differently colored.
#' @param open.angle Angle aperture circular layout. Default = \code{100}.
#' @param scale.size Numerical value indicating how much scale points in the dendogram terminal nodes. Default = \code{10}.
#' @param colors Vector with colors for plotting.
#' @param labels Two values indicating if labels should be displayed (default, \code{TRUE}) and the offset value (for separating labels from tree tips, 0.2 as default).
#' @keywords circular tree
#' @keywords dendrogram
#' @keywords nodes
#' @keywords hierachical clustering
#' @export
#' @import dplyr
#' @import ggplot2
#' @import ggtree
#' @importFrom stats dist hclust median
#' @examples
#' \dontrun{
#' # step 1: display all node numbers to select how to coloring areas
#' circ.tree(fcs.SCE = fcsL, cell.clusters = "SOM_L_named", nodes = "display")
#' 
#' # step 2: color areas indicating node numbers
#' circ.tree(fcs.SCE = fcsL, cell.clusters = "SOM_L_named", nodes = c(17, 25))
#' }

circ.tree <- function(fcs.SCE, assay.i = "normalized", cell.clusters, dist.method = "euclidean", hclust.method = "average", 
                      nodes = "display", open.angle = 50, scale.size = 10, colors, labels = c(T, 0.2)){
  if (!requireNamespace(c("ape", "ggtree"), quietly = TRUE)){
    cat("Package \"ggtree\" is needed for this function, installing it...\n\n")
    BiocManager::install(c("ape", "ggtree")); require(ggtree)
  } else require(ggtree)
  
  exprs <- t(SummarizedExperiment::assay(fcs.SCE, i = assay.i))
  exprs_01 <- scale.exprs(exprs)
  if(missing(colors)) colors <- div.colors(length(unique(fcs.SCE[[cell.clusters]])))
  
  ## median expression
  expr_medianL <- data.frame(exprs, cell_clustering = fcs.SCE[[cell.clusters]]) %>%
    group_by(.data$cell_clustering) %>% dplyr::summarize_all(list(median)) %>% as.data.frame(.data)
  colnames(expr_medianL)[1] <- cell.clusters
  
  expr01_medianL <- data.frame(exprs_01, cell_clustering = fcs.SCE[[cell.clusters]]) %>%
    group_by(.data$cell_clustering) %>% dplyr::summarize_all(list(median)) %>% as.data.frame(.data)
  colnames(expr01_medianL)[1] <- cell.clusters
  
  ## calculate cluster frequencies
  clustering_propL <- data.frame(node = factor(1:length(unique(fcs.SCE[[cell.clusters]]))), prop.table(table(fcs.SCE[[cell.clusters]]))*100)
  colnames(clustering_propL) <- c("label", "cell_cluster", "Freq")
  
  ## hierarchical clustering
  suppressWarnings(dL <- dist(expr_medianL, method = dist.method))
  cluster_rowsL <- hclust(dL, method = hclust.method)
  
  hca <- ape::as.phylo(cluster_rowsL)
  # hca <- full_join(hca, data.frame(label = as.character(1:nrow(expr01_medianL)), 
  #                                  label2 = expr01_medianL[,cell.clusters]), by = "label")
  hca <- dplyr::full_join(hca, clustering_propL, by = "label")
  
  ## plotting
  if(length(nodes) == 1 && nodes == "display"){
    return(ggtree::ggtree(hca, layout = "fan") + ggtree::geom_text2(aes(label = node), hjust = -.3) + ggtree::geom_tiplab(aes(label = cell_cluster)))
  }else{
    p1 <- ggtree::ggtree(hca, layout = "fan", open.angle = open.angle, root.position = 10)
    
    for(i in nodes){
      p1 <- p1 + ggtree::geom_hilight(node = i, fill = sample(colors, 1), alpha = .6) +
        geom_point(aes_string(color = "cell_cluster", size = "Freq")) +
        scale_size_area(max_size = scale.size) + scale_color_manual(values = colors) + 
        theme(legend.position = "bottom")
    }
  }
  
  g1 <- ggtree::gheatmap(p1, expr01_medianL[,-1], offset = 0.01, width = 1, font.size = 4, colnames_angle = 45, hjust = 0,
                         colnames_position = "top", high = "#b30000", low = "#fff7f3") + theme(legend.position = "bottom")
  if(labels[1]) return(g1 + geom_tiplab(aes(label = cell_cluster), offset = labels[2], align = T)) else return(g1)
}


## from gtree source...
`%add%` <- function(p, data) {
  p$data <- p$data %add2% data
  return(p)
}

`%add2%` <- function(d1, d2) {
  if ("node" %in% colnames(d2)) {
    cn <- colnames(d2)
    ii <- which(cn %in% c("node", cn[!cn %in% colnames(d1)]))
    d2 <- d2[, ii]
    dd <- dplyr::left_join(d1, d2, by="node")
  } else {
    d2[,1] <- as.character(unlist(d2[,1])) ## `unlist` to work with tbl_df
    d2 <- dplyr::rename(d2, label = 1) ## rename first column name to 'label'
    dd <- dplyr::left_join(d1, d2, by="label")
  }
  dd <- dd[match(d1$node, dd$node),]
  return(dd)
}
