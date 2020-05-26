# 'sc.Heatmap
#'
#' This function draws a heatmap with single-cell fluorescence values.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param not.metadata Vector with variable names (from \code{colData(fcs.SCE)}) for not including in the heatmap annotation. Default = \code{c("filename", "cell_id", "sample_id")}.
#' @param clustering.method Clustering method for rows and columns clustering within the heatmap. Possible values are "average" (default), "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid".
#' @param subsampling Numeric value indicating how many events use to draw heatmap and speed up plotting. Default = \code{100}.
#' @param color Color vector for coloring expression values within the heatmap. Default = \code{NULL} (i.e., scale "YlGnBu" from \code{RColorBrewer}).
#' @keywords single-cell expression
#' @keywords heatmap
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom SummarizedExperiment colData assay
#' @importFrom pheatmap pheatmap
#' @examples
#' \dontrun{
#' sc.heatmap(fcs, subsampling = 100)
#' }

sc.heatmap <- function(fcs.SCE, assay.i = "normalized", markers.to.use = "all", not.metadata = c("filename", "cell_id", "sample_id"), 
                       clustering.method = "average", subsampling = 100, color = NULL){
  if(!is.null(subsampling)) suppressMessages(fcs.SCE <- sub.samples(fcs.SCE, subsampling = subsampling))
  if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- rownames(fcs.SCE)
  if(is.null(color)) color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
  data <- assay(fcs.SCE, i = assay.i)
  metadata <- as.data.frame(colData(fcs.SCE))
    
  annot_col <- col.annot.pheatmap(metadata[,!(colnames(metadata) %in% not.metadata)])

  pheatmap(data[markers.to.use2,], 
           annotation_col = metadata[,!(colnames(metadata) %in% not.metadata)], 
           annotation_colors = annot_col, 
           color = color, display_numbers = FALSE, 
           number_color = "black", fontsize_number = 5, 
           show_colnames = F, clustering_method = clustering.method, 
           treeheight_row = 0, treeheight_col = 0)
}