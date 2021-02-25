#' Heatmap from median expression values
#'
#' This function draws a heatmap with median values for each FCS file or for identified cluster with \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}}
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters Name of column containing clusters identified through \code{\link[FlowCT:clustering.flow]{FlowCT::clustering.flow()}}.
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param not.metadata Vector with variable names (from \code{colData(fcs.SCE)}) for not including in the heatmap annotation. Default = \code{"filename"}.
#' @param colors Vector with colors for plotting (if provided, it must be as long as the number of unique elements in the longer metadata field). Default = \code{NULL} (i.e., it will choose automatically a vector of colors according to \code{\link[FlowCT:div.colors]{FlowCT::div.colors()}}).
#' @keywords heatmap
#' @keywords cell cluster percentages
#' @keywords median expression values
#' @export median.heatmap
#' @import dplyr
#' @importFrom stats median dist hclust
#' @importFrom SummarizedExperiment assay
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' # option 1: general heatmap (by FCS file)
#' median.heatmap(fcs, not.metadata = c("sample_id", "file_name"))
#' 
#' # option 2: heatmap by SOM-identified clusters
#' median.heatmap(fcs.SCE = fcs, assay.i = "normalized", cell.clusters = "SOM")
#' }

median.heatmap <- function(fcs.SCE, assay.i = "normalized", cell.clusters = NULL, markers.to.use = "all", not.metadata = "filename", colors = NULL){
  data <- t(assay(fcs.SCE, i = assay.i))
  colnames(data) <- make.names(colnames(data)) #avoid R renaming conflicts
  metadata <- dplyr::distinct(as.data.frame(colData(fcs.SCE)), .data$filename, .keep_all = T)
  
  if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- make.names(colnames(data))
  
  ## prepare median tables
  if(is.null(cell.clusters)){
    med <- median.values(fcs.SCE, assay.i = assay.i)
  }else{
    expr_median <- data.frame(cell_clusters = metadata[,cell.clusters], data[,markers.to.use]) %>%
      group_by(.data$cell_clusters) %>% summarize_all(list(median)) %>% as.data.frame(.data)
    
    expr_saturated_median <- data.frame(cell_clusters = metadata[,cell.clusters], scale.exprs(data[,markers.to.use])) %>%
      group_by(.data$cell_clusters) %>% summarize_all(list(median)) %>% as.data.frame(.data)
  }
  
  ## heatmap
  if(is.null(cell.clusters)){
    annotation_colors <- col.annot.pheatmap(metadata[,!(colnames(metadata) %in% not.metadata), drop = F], colors = colors)
    color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
    
    print(pheatmap::pheatmap(t(med[,markers.to.use]), color = color, display_numbers = FALSE,
                   number_color = "black", fontsize_number = 5, clustering_method = "average",
                   annotation = metadata[,!(colnames(metadata) %in% not.metadata), drop = F], 
                   annotation_colors = annotation_colors, 
                   show_colnames = F))
  }else{
    ## calculate cluster frequencies
    clustering_table <- table(cell.clusters)
    clustering_prop <- round(prop.table(table(metadata[,cell.clusters]))*100, digits = 2)
    labels_row <- paste0(expr_saturated_median$cell_clusters, " (", clustering_prop ,"%)")
    
    d <- dist(expr_median[,markers.to.use], method = "euclidean")
    cluster_rows <- hclust(d, method = "average")
    expr_heat <- as.matrix(expr_saturated_median[,markers.to.use])
    # rownames(expr_heat) <- paste0("c.", expr_saturated_median$cell_clusters) #force rownames to not crash heatmap (??)
    rownames(expr_heat) <- rownames(expr_saturated_median) #force rownames to not crash heatmap (??)
    
    ## annot colors
    annot_row <- expr_saturated_median[,"cell_clusters", drop = F]
    # rownames(annot_row) <- paste0("c.", rownames(annot_row)) #force rownames to not crash heatmap (??)
    annotation_colors <- col.annot.pheatmap(expr_saturated_median[,"cell_clusters", drop = F], colors = colors)
    color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
    legend_breaks <- seq(from = 0, to = 1, by = 0.1)
    
    print(pheatmap::pheatmap(expr_heat, color = color, annotation_legend = F,
                   cluster_cols = FALSE, cluster_rows = cluster_rows, labels_row = labels_row,
                   display_numbers = FALSE, number_color = "black",
                   fontsize = 8, fontsize_number = 6, legend_breaks = legend_breaks,
                   annotation_row = annot_row, annotation_colors = annotation_colors))
  }
}
