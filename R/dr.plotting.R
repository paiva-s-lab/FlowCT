#' dr.plotting
#'
#' This function plots the indicated dimensional reduction (DR) from a previously calculated \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}} object.
#' @param dr A object with DR generated with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}} or a \code{data.frame} with DR, expression and metadata information (like the first element list of the object generated with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}}).
#' @param dr.calculated String indicating the desired DR to plot (this indicated DR should be prevoulsy calculated to being plotted).
#' @param n.dims Vector indicating the two DR components to plot. Default = \code{c(1,2)}.
#' @param color.by Variable from (from \code{colData(fcs.SE)}) for dots coloring. If \code{color.by = "expression"} (default), plot will be splitted for each marker (\code{facet.by}).
#' @param facet.by Variable from (from \code{colData(fcs.SE)}) for plot spliting. Default = \code{NULL}.
#' @param omit.markers Vector with markers to omit when plotting with \code{color.by = "expression"}. Default = \code{NULL}.
#' @param title Title to add to the plot.
#' @param label.by Variable from (from \code{colData(fcs.SE)}) for dots labeling. Default = \code{NULL}.
#' @param size Point size. Default = \code{0.5}.
#' @param raster Vector indicating if image should be rasterized (logical element) and the number of pixels to consider (numerical element). It is based on (\ref{https://github.com/exaexa/scattermore}{scattermore package}). Default = \code{c(T, 1000)}.
#' @keywords dimensional reduction plotting
#' @keywords tSNE
#' @keywords PCA
#' @keywords UMAP
#' @export
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' dr.plotting(dr, dr.calculated = "tSNE", color.by = "condition")
#' dr.plotting(dr, dr.calculated = "UMAP", color.by = "patient_id")
#' dr.plotting(dr, dr.calculated = "PCA", color.by = "SOM", facet.by = "condition")
#' }

dr.plotting <- function(dr, dr.calculated, n.dims = c(1,2), color.by = "expression", facet.by = NULL, omit.markers = NULL, title = "", label.by = NULL, size = 2, raster = c(T, 1000)){
  require(ggplot2)

  if(class(dr) == "list"){
    if(is.null(omit.markers)){
      if(color.by == "expression"){
        dr <- dr$dr_melted
      }else{
        dr <- dr$dr
      }
    }else{
      if(color.by == "expression"){
        dr <- dr$dr_melted[!grepl(paste0(omit.markers, collapse = "|"), dr$dr_melted$antigen),]
      }else{
        dr <- dr$dr[,!grepl(paste0(omit.markers, collapse = "|"), colnames(dr$dr))]
      }
    }
  }else{
    dr <- dr
  }
  
  g <- ggplot(dr, aes_string(x = paste0(dr.calculated, n.dims[1]), y = paste0(dr.calculated, n.dims[2]), 
                             color = color.by)) +
    xlab(paste0(dr.calculated, 1)) + ylab(paste0(dr.calculated, 2)) +
    ggtitle(title) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA))
  
  if(is.factor(dr[,color.by]) | is.character(dr[,color.by])){
    g <- g + scale_color_manual(values = div.colors(length(unique(dr[,color.by]))), name = color.by) +
      guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
  }else{
    g <- g + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50), name = color.by)
  }
  
  if(!is.null(facet.by)){
    g <- g + facet_wrap(~ eval(parse(text = facet.by)))
  }else if(color.by == "expression" & is.null(facet.by)){
    g <- g + facet_wrap(~ antigen)
  }

  if(!is.null(label.by)){
    g <- g + geom_text(aes_string(label = label.by), nudge_y = 0.05)
  }

  if(raster[1]){
    g <- g + scattermore::geom_scattermore(pointsize = size, pixels = rep(raster[2], 2)) #devtools::install_github('exaexa/scattermore')
  }else{
    g <- g + geom_point(aes_string(color = color.by), size = size)
  }
  
  return(g)
}
