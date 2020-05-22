#' cell.count.bx
#'
#' This function draws a boxplot according to the cell count for a specified condition.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param x.axis Variable name (from \code{colData(fcs.SCE)}) for drawing the x-axis in plot.
#' @param color.by Variable name (from \code{colData(fcs.SCE)}) for coloring the plot. Default = \code{x.axis} (i.e., the same specfied in \code{x.axis}).
#' @param label.by Variable from (from \code{colData(fcs.SCE)}) for labeling. Default = \code{NULL}.
#' @param limits Numeric vector with limits for plotting (minimum, maximum). Default = \code{NULL} (i.e., automatically calculated from data).
#' @keywords cell count
#' @keywords boxplot
#' @export
#' @import ggplot2
#' @importFrom SummarizedExperiment colData assay
#' @importFrom ggpubr stat_compare_means ggboxplot
#' @importFrom ggrepel geom_label_repel
#' @examples
#' \dontrun{cell.count.bx(fcsL, assay.i = "normalized", x.axis = "condition")}

cell.count.bx <- function(fcs.SCE, assay.i = "normalized", x.axis, color.by = x.axis, label.by = NULL, limits = NULL){
  data <- merge(t(assay(fcs.SCE, i = assay.i)), colData(fcs.SCE), by = "row.names")[,-1]
  ggdf <- data.frame(data[!duplicated(data[,"filename"]),], cell_counts = as.numeric(table(data[,"filename"])))
  
  if(is.null(limits)) limits <- c(min(ggdf$cell_counts), max(ggdf$cell_counts))
  
  g <- ggboxplot(ggdf, x = x.axis, y = "cell_counts", fill = color.by) +
    stat_compare_means(label.x = 1.7, label.y = max(ggdf$cell_counts)) +
    scale_y_continuous(limits = limits) + 
    scale_fill_manual(values = div.colors(length(unique(ggdf[,color.by]))), drop = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    geom_jitter(shape=16, position=position_jitter(0.2))
  
  if(!is.null(label.by)){
    print(g + geom_label_repel(aes(label = eval(parse(text = label.by))), alpha=0.75, fontface = 'bold', color = 'black'))
  }else{
    print(g)
  }
}