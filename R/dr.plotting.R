#' Plot DR data
#'
#' This function plots the indicated dimensional reduction (DR) from a previously calculated \code{\link[FlowCT.v2:dim.reduction]{FlowCT.v2::dim.reduction()}} object.
#' @param data A object with DR generated with \code{\link[FlowCT.v2:dim.reduction]{FlowCT.v2::dim.reduction()}} or a \code{data.frame} with DR, expression and metadata information (like the first element list of the object generated with \code{\link[FlowCT.v2:dim.reduction]{FlowCT.v2::dim.reduction()}}).
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param plot.dr String indicating the desired DR to plot (this indicated DR should be prevoulsy calculated to being plotted).
#' @param n.dims Vector indicating the two DR components to plot. Default = \code{c(1,2)} (by now, these are the only dims allowed).
#' @param color.by Variable from (from \code{colData(fcs.SCE)}) for dots coloring. If \code{color.by = "expression"} (default), plot will be splitted for each marker (\code{facet.by}).
#' @param shape.by Variable from (from \code{colData(fcs.SCE)}) for dots shaping. Default = \code{NULL}.
#' @param facet.by Variable from (from \code{colData(fcs.SCE)}) for plot spliting. Default = \code{NULL}.
#' @param omit.markers Vector with markers to omit when plotting with \code{color.by = "expression"}. Default = \code{NULL}.
#' @param title Title to add to the plot.
#' @param label.by Variable from (from \code{colData(fcs.SCE)}) for dots labeling. Default = \code{NULL}.
#' @param size Point size. Default = \code{1}.
#' @param raster Vector indicating if image should be rasterized (logical element) and the number of pixels to consider (numerical element). It is based on \href{https://github.com/exaexa/scattermore}{\code{scattermore} package}. Default = \code{c(FALSE, 1000)}.
#' @param return.df Logical indicating if built \code{data.frame} with DR information and metadata must be returned. Default = \code{FALSE}.
#' @param colors Vector with colors for plotting. Default = \code{NULL} (i.e., it will choose automatically a vector of colors according to \code{\link[FlowCT.v2:div.colors]{FlowCT.v2::div.colors()}}).
#' @keywords dimensional reduction plotting
#' @keywords tSNE
#' @keywords PCA
#' @keywords UMAP
#' @export
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom data.table melt as.data.table
#' @importFrom scattermore geom_scattermore
#' @examples
#' \dontrun{
#' dr.plotting(fcs, plot.dr = "tSNE", color.by = "condition")
#' dr.plotting(fcs, plot.dr = "UMAP", color.by = "patient_id")
#' dr <- dr.plotting(fcs, plot.dr = "PCA", color.by = "SOM", facet.by = "condition", return.df = T)
#' }

dr.plotting <- function(data, assay.i = "normalized", plot.dr, dims = c(1,2), color.by = "expression", shape.by = NULL, facet.by = NULL, omit.markers = NULL, title = "", label.by = NULL, size = 1, raster = c(F, 1000), return.df = F, colors = NULL){
  if(class(data)[1] == "SingleCellExperiment"){
    pos <- match(tolower(plot.dr), tolower(names(data@int_colData@listData$reducedDims)))
    if(is.na(pos)) stop('The DR indicated has not been calculated yet or is differently named (please, check the output of reducedDimNames(data) to see the correct DR name).\n', call. = F)
    dr_calculated <- names(data@int_colData@listData$reducedDims)[pos]
    dr <- data@int_colData@listData$reducedDims@listData[[dr_calculated]][,dims]
    colnames(dr) <- paste0("dr", dims) #useless

    no.omit.markers <- rownames(data)[!(rownames(data) %in% omit.markers)]
    drmd <- as.data.frame(cbind(colData(data), dr, t(assay(data, i = assay.i))[,no.omit.markers]))
  }else{
    drmd <- data
  }

  if(color.by == "expression") drmd <- as.data.frame(melt(as.data.table(drmd), measure.vars = no.omit.markers, value.name = "expression", variable.name = "antigen"))
  if(is.null(colors)) colors <- div.colors(length(unique(drmd[,color.by])))

  g <- ggplot(drmd, aes_string(x = paste0("dr", dims[1]), y = paste0("dr", dims[2]), color = color.by)) +
    xlab(paste0(toupper(plot.dr), "-1")) + ylab(paste0(toupper(plot.dr), "-2")) + ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA))

  if(raster[1]){
    g <- g + geom_scattermore(pointsize = size, pixels = rep(raster[2], 2)) #devtools::install_github('exaexa/scattermore')
  }else{
    g <- g + geom_point(aes_string(color = color.by), size = size)
  }

  if(is.numeric(drmd[,color.by])){
    g <- g + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50), name = color.by)
  }else{
    g <- g + scale_color_manual(values = colors, name = color.by) +
      guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
  }

  if(!is.null(facet.by)){
    g <- g + facet_wrap(~ eval(parse(text = facet.by)))
  }else if(color.by == "expression" & is.null(facet.by)){
    g <- g + facet_wrap(~ antigen)
  }

  if(!is.null(label.by)){
    g <- g + geom_text(aes_string(label = label.by), nudge_y = 0.05)
  }

  print(g)
  if(return.df) return(drmd)
}

