#' median.heatmap
#'
#' This function draws a heatmap with median values for each FCS file or for identified cluster with \code{\link[FlowCT.v2:fsom.clustering]{FlowCT.v2::fsom.clustering()}}
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT.v2:fsom.clustering]{FlowCT.v2::fsom.clustering()}} (and, normaly, later renamed). Default = \code{NULL} (i.e., median values will be calculated for each FCS file).
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param not.metadata Vector with variable names (from \code{colData(fcs.SCE)}) for not including in the heatmap annotation. Default = \code{"filename"}.
#' @keywords heatmap
#' @keywords cell cluster percentages
#' @keywords median expression values
#' @export median.heatmap
#' @import dplyr
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats median dist hclust
#' @importFrom pheatmap pheatmap
#' @examples
#' \dontrun{
#' # option 1: general heatmap (by FCS file)
#' median.heatmap(fcs, not.metadata = c("sample_id", "file_name"))
#' 
#' # option 2: heatmap by SOM-identified clusters
#' median.heatmap(fcs.SCE = fcs, assay.i = "normalized", cell.clusters = fcs$SOM)
#' }

median.heatmap <- function(fcs.SCE, assay.i = "normalized", cell.clusters = NULL, markers.to.use = "all", not.metadata = "filename"){
  data <- t(assay(fcs.SCE, i = assay.i))
  metadata <- fcs.SCE@metadata$reduced_metadata
  
  if(markers.to.use == "all") markers.to.use2 <- colnames(data) else markers.to.use2 <- markers.to.use
  
  ## prepare median tables
  if(is.null(cell.clusters)){
    med <- median.values(fcs.SCE, assay.i = assay.i)
  }else{
    expr_median <- data.frame(cell_clusters = cell.clusters, data[,markers.to.use2]) %>%
      group_by(.data$cell_clusters) %>% summarize_all(list(median)) %>% as.data.frame(.data)
    
    expr_saturated_median <- data.frame(cell_clusters = cell.clusters, scale.exprs(data[,markers.to.use2])) %>%
      group_by(.data$cell_clusters) %>% summarize_all(list(median)) %>% as.data.frame(.data)
  }
  
  ## heatmap
  if(is.null(cell.clusters)){
    annotation_colors <- col.annot.pheatmap(metadata[,!(colnames(metadata) %in% not.metadata), drop = F])
    color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
    
    print(pheatmap(t(med[,markers.to.use2]), color = color, display_numbers = FALSE,
                   number_color = "black", fontsize_number = 5, clustering_method = "average",
                   annotation = metadata[,!(colnames(metadata) %in% not.metadata), drop = F], 
                   annotation_colors = annotation_colors, 
                   show_colnames = F))
  }else{
    ## calculate cluster frequencies
    clustering_table <- table(cell.clusters)
    clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
    labels_row <- paste0(expr_saturated_median$cell_clusters, " (", clustering_prop ,"%)")
    
    d <- dist(expr_median[,markers.to.use2], method = "euclidean")
    cluster_rows <- hclust(d, method = "average")
    expr_heat <- as.matrix(expr_saturated_median[,markers.to.use2])
    # rownames(expr_heat) <- paste0("c.", expr_saturated_median$cell_clusters) #force rownames to not crash heatmap (??)
    rownames(expr_heat) <- rownames(expr_saturated_median) #force rownames to not crash heatmap (??)
    
    ## annot colors
    annot_row <- expr_saturated_median[,"cell_clusters", drop = F]
    # rownames(annot_row) <- paste0("c.", rownames(annot_row)) #force rownames to not crash heatmap (??)
    annotation_colors <- col.annot.pheatmap(expr_saturated_median[,"cell_clusters", drop = F])
    color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
    legend_breaks <- seq(from = 0, to = 1, by = 0.1)
    
    print(pheatmap(expr_heat, color = color, annotation_legend = F,
                   cluster_cols = FALSE, cluster_rows = cluster_rows, labels_row = labels_row,
                   display_numbers = FALSE, number_color = "black",
                   fontsize = 8, fontsize_number = 6, legend_breaks = legend_breaks,
                   annotation_row = annot_row, annotation_colors = annotation_colors))
  }
}
