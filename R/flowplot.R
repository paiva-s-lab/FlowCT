#' Draw expression dotplots
#'
#' It draws the classical flow cytometry dotplot faceting two markers. Multiple dot plots can be plotted if \code{x.axis} and \code{y.axis} are specified with multiple markers.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}. By default, the matrix used is the \code{arcsinh} transformed.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param x.axis Vector of markers to draw x-axis on the dotplots. Length must be the same of \code{y.axis} because each marker has its parter in \code{y.axis}.
#' @param y.axis Vector of markers to draw y-axis on the dotplots. Length must be the same of \code{x.axis} because each marker has its parter in \code{x.axis}.
#' @param densities Logical indicating if densities must be drawn in each plot. Default = \code{TRUE}.
#' @param color.by Variable from (from \code{colData(fcs.SCE)}) for events coloring.
#' @param select.values.color Vector of values taken from \code{color.by} option to include in the coloring. Default = \code{"all"}.
#' @param size Point (event) size. Default = \code{0.5}.
#' @param colors Vector with colors for plotting. Default = \code{NULL} (i.e., it will choose automatically a vector of colors according to \code{\link[FlowCT.v2:div.colors]{FlowCT.v2::div.colors()}}).
#' @keywords scatterplot
#' @keywords dotplot
#' @keywords marker
#' @export
#' @import ggplot2
#' @import dplyr
#' @import cowplot
#' @importFrom gridExtra grid.arrange
#' @importFrom SummarizedExperiment assay colData
#' @examples
#' \dontrun{
#' flowplot(fcs, x.axis = c("FSC_A", "CD4", "CD62L"), y.axis = c("SSC_A", "CD8", "CD45"),
#'     color.by = "SOM", select.values.color = 1:10)
#' }

flowplot <- function(fcs.SCE, assay.i = "normalized", x.axis, y.axis, densities = T, color.by, select.values.color = "all", size = 0.5, colors = NULL){
  data <- as.data.frame(cbind(colData(fcs.SCE), t(assay(fcs.SCE, assay.i))))

  ## filter data
  if(length(select.values.color) == 1 && select.values.color == "all") select.values.color <- unique(data[,color.by])
  data <- data[data[,color.by] %in% select.values.color,]

  ## plotting
  if(is.null(colors)) colors <- div.colors(length(select.values.color))
  combos <- cbind(as.data.frame(x.axis), as.data.frame(y.axis))

  gglist <- list()
  for(i in 1:nrow(combos)){
    g1 <- ggplot(data, aes_string(x = combos[i,1], y = combos[i,2], color = color.by)) +
      geom_point(size = size) +
      scale_color_manual(values = colors) +
      theme_bw() + guides(color = guide_legend(ncol = 3))
    legend <- get_legend(g1)

    if(densities){
      xdens <- axis_canvas(g1, axis = "x")+ #histogram along x axis
        scale_fill_manual(values = colors) +
        geom_density(data = data, aes_string(x = combos[i,1], fill = color.by),
                     alpha = 0.7, size = 0.2)

      ydens <- axis_canvas(g1, axis = "y", coord_flip = TRUE) + #histogram along y axis
        scale_fill_manual(values = colors) +
        geom_density(data = data, aes_string(x = combos[i,2], fill = color.by), alpha = 0.7, size = 0.2)+
        coord_flip()

      p1 <- insert_xaxis_grob(g1 + theme(legend.position = "none"), xdens, grid::unit(.2, "null"), position = "top")
      gglist[[i]] <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
    }else{
      gglist[[i]] <- g1 + theme(legend.position = "none")
    }
  }
  gglist$legend <- legend
  do.call(gridExtra::grid.arrange, gglist)
}
