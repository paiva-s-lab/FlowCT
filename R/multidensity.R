#' multidensity
#'
#' This function draws multidensity plot with all FCS files included in a \code{FCS.SE} object. If there are more files than limit specified in \code{ridgeline.lim}, instead of plottting density in a ridge-way all density lines will be overlapped.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param show.markers Vector with markers to plot. Default = \code{"all"}.
#' @param color.by Variable name (from \code{colData(fcs.SE)}) for lines coloring. 
#' @param subsampling Numeric value indicating how many events use to calculate density lines and speed up plotting. Default = \code{NULL}.
#' @param interactive Logical indicating if the user can interact with the (only overlapping-lines) plot. Default = \code{FALSE}.
#' @param ridgeline.lim Numeric value specifying the limit for shifting from ridgeline-mode to overlapping-lines plot. Default = \code{15}.
#' @keywords marker normalization
#' @keywords marker density
#' @keywords marker alignment
#' @export
#' @import ggplot2
#' @importFrom SummarizedExperiment colData assay
#' @importFrom stats median
#' @examples
#' \dontrun{
#' multidensity(fcs.SE = fcs_se, assay.i = "normalized", subsampling = 1000)
#' multidensity(fcs_se, assay.i = 2, color.by = "file_name", ridgeline.lim = 0, show.markers = c("CD62L", "CD4"), interactive = T)
#' }

multidensity <- function(fcs.SE, assay.i, show.markers = "all", color.by = NULL, subsampling = NULL, interactive = F, ridgeline.lim = 15){
  if(show.markers == "all") show.markers <- rownames(fcs.SE)
  if(!is.null(subsampling)) suppressMessages(fcs.SE <- sub.samples(fcs.SE, subsampling = subsampling))
  
  data <- t(assay(fcs.SE, i = assay.i))
  data2 <- merge(data, colData(fcs.SE), by = "row.names")[,-1]
  
  # prepare tables: for plotting and with median values for each marker
  median_df <- data.frame(antigen = show.markers, median = apply(data[,show.markers], 2, median))
  ggdf <- data.table::melt(data2, measure.vars = show.markers, value.name = "expression", variable.name = "antigen")
  
  if(length(unique(fcs.SE$filename)) > ridgeline.lim){
    g <- ggplot(data = ggdf[grepl(paste0(show.markers, collapse = "|"), ggdf$antigen),], 
                aes_string(x = "expression", color = color.by, group = "filename")) + 
      # geom_density(size = 0.5) +
      stat_density(geom = "line", position = "identity", size = 0.5) +
      facet_wrap(~ antigen, scales = "free") +
      geom_vline(data = median_df, aes(xintercept = median), linetype = 2, color = "gray55") +
      scale_color_manual(name = color.by, values = div.colors(unique(length(ggdf[,color.by])))) +
      theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                              strip.text = element_text(size = 7), axis.text = element_text(size = 5))
    
    if(interactive) plotly::ggplotly(g) else print(g)
  }else{
    suppressMessages(print(ggplot(ggdf, aes_string(x = "expression", y = "filename")) + 
                             ggridges::geom_density_ridges(alpha = 0.7) +
                             facet_wrap(~ antigen, scales = "free") +
                             geom_vline(data = median_df, aes(xintercept = median), linetype = 2, color = "gray55") +
                             theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                                     strip.text = element_text(size = 7), axis.text = element_text(size = 7))))
  }
}