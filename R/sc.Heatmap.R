# 'sc.Heatmap
#'
#' This function draws a heatmap with single-cell fluorescence values.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param not.metadata Vector with variable names (from \code{colData(fcs.SE)}) for not including in the heatmap annotation. Default = \code{c("filename", "cell_id", "sample_id")}.
#' @param clustering.method Clustering method for rows and columns clustering within the heatmap. Possible values are "average" (default), "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid".
#' @param subsampling Numeric value indicating how many events use to draw heatmap and speed up plotting. Default = \code{100}.
#' @param color Color vector for coloring expression values within the heatmap. Default = \code{NULL} (i.e., scale "YlGnBu" from \code{RColorBrewer}).
#' @keywords single-cell expression
#' @keywords heatmap
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom SummarizedExperiment colData assay
#' @examples
#' \dontrun{
#' sc.heatmap(fcs_se100)
#' }

sc.heatmap <- function(fcs.SE, assay.i = "normalized", markers.to.use = "all", not.metadata = c("filename", "cell_id", "sample_id"), 
                       clustering.method = "average", subsampling = 100, color = NULL){
  if(!is.null(subsampling)) suppressMessages(fcs.SE <- sub.samples(fcs.SE, subsampling = subsampling))
  if(markers.to.use == "all") markers.to.use <- rownames(fcs.SE)
  data <- t(assay(fcs.SE, i = assay.i))
  metadata <- as.data.frame(colData(fcs.SE))
    
  annot_col <- col.annot.pheatmap(metadata[,!(colnames(metadata) %in% not.metadata)])
  
  pheatmap::pheatmap(t(data[,markers.to.use]), 
           annotation_col = metadata[,!(colnames(metadata) %in% not.metadata)], 
           annotation_colors = annot_col, color = color, display_numbers = FALSE, 
           number_color = "black", fontsize_number = 5, 
           show_colnames = F, clustering_method = clustering.method, 
           treeheight_row = 0, treeheight_col = 0)

}