#' Density plots with expression values overlapped
#'
#' This function draws multidensity plot with all FCS files included in a \code{fcs.SCE} object. If there are more files than limit specified in \code{ridgeline.lim}, instead of plottting density in a ridge-way all density lines will be overlapped.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param show.markers Vector with markers to plot. Default = \code{"all"}.
#' @param color.by Variable name (from \code{colData(fcs.SCE)}) for lines coloring.
#' @param subsampling Numeric value indicating how many events use to calculate density lines and speed up plotting. Default = \code{NULL}.
#' @param interactive Logical indicating if the user can interact with the (only overlapping-lines) plot. Default = \code{FALSE}.
#' @param ridgeline.lim Numeric value specifying the limit for shifting from ridgeline-mode to overlapping-lines plot. Default = \code{15}.
#' @param colors Vector with colors for plotting (only available for no-ridgeline plot). Default = \code{NULL} (i.e., it will choose automatically a vector of colors according to \code{\link[FlowCT:div.colors]{FlowCT::div.colors()}}).
#' @keywords marker normalization
#' @keywords marker density
#' @keywords marker alignment
#' @export
#' @import ggplot2
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom stats median
#' @importFrom data.table melt as.data.table
#' @examples
#' \dontrun{
#' multidensity(fcs.SCE = fcs, assay.i = "normalized", subsampling = 1000)
#' multidensity(fcs, assay.i = 2, color.by = "file_name", ridgeline.lim = 0, 
#'      show.markers = c("CD62L", "CD4"), interactive = T)
#' }

multidensity <- function(fcs.SCE, assay.i, show.markers = "all", color.by = NULL, subsampling = NULL, interactive = F, ridgeline.lim = 15, colors = NULL){
  if(show.markers == "all") show.markers <- make.names(rownames(fcs.SCE)) else show.markers <- make.names(show.markers)
  if(!is.null(subsampling)) suppressMessages(fcs.SCE <- sub.samples(fcs.SCE, subsampling = subsampling))

  data <- t(assay(fcs.SCE, i = assay.i))
  colnames(data) <- make.names(colnames(data))
  data2 <- cbind(data, colData(fcs.SCE))

  # prepare tables: for plotting and with median values for each marker
  median_df <- data.frame(antigen = show.markers, median = apply(data[,show.markers], 2, median))
  ggdf <- as.data.frame(data.table::melt(data.table::as.data.table(data2), measure.vars = show.markers, value.name = "expression", variable.name = "antigen"))

  if(is.null(colors)) colors <- div.colors(unique(length(ggdf[,color.by])))

  if(length(unique(fcs.SCE$filename)) > ridgeline.lim){
    g <- ggplot(data = ggdf[grepl(paste0(show.markers, collapse = "|"), ggdf$antigen),],
                aes_string(x = "expression", color = color.by, group = "filename")) +
      # geom_density(size = 0.5) +
      stat_density(geom = "line", position = "identity", size = 0.5) +
      facet_wrap(~ antigen, scales = "free") +
      geom_vline(data = median_df, aes(xintercept = median), linetype = 2, color = "gray55") +
      scale_color_manual(name = color.by, values = colors) +
      theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                              strip.text = element_text(size = 7), axis.text = element_text(size = 5))

    if(interactive){
      if(!requireNamespace("plotly", quietly = TRUE)) stop("Package \"plotly\" needed for this function to work. Please install it.", call. = FALSE)
      plotly::ggplotly(g)
    }else return(g)
  }else{
    if(!requireNamespace("ggridges", quietly = TRUE)) stop("Package \"ggridges\" needed for this function to work. Please install it.", call. = FALSE)
    return(ggplot(ggdf, aes_string(x = "expression", y = "filename", fill = color.by)) +
             ggridges::geom_density_ridges(alpha = 0.7, color = "gray25") +
             facet_wrap(~ antigen, scales = "free") +
             geom_vline(data = median_df, aes(xintercept = median), linetype = 2, color = "gray55") +
             theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                     strip.text = element_text(size = 7), axis.text = element_text(size = 7)))
  }
}
