#' Boxplots for identified clusters
#'
#' It draws a boxplot with cell clusters identified through \code{\link[FlowCT.v2:fsom.clustering]{FlowCT.v2::fsom.clustering()}}.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT.v2:fsom.clustering]{FlowCT.v2::fsom.clustering()}} (and, normaly, later renamed).
#' @param condition.column Column name from the \code{colData(fcs.SCE)} object which contains condition information.
#' @param pvalue.cutoffs List of P-value cutoffs and their symbols for indicante significances within the plot. Default = \code{list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("xxxx", "***", "**", "*", "ns"))}.
#' @param color.by Variable name (from \code{colData(fcs.SCE)}) for lines coloring. 
#' @param geom.point Logical indicating if adding points to boxplot. Default = \code{TRUE}.
#' @param facet Logical indicating if splitting boxplots by cell clusters. Default = \code{FALSE}.
#' @param facet.free.scale If \code{facet = TRUE}, string indicating how scales would be shared across all facets. Possible values: \code{"free_x"} (default), \code{"free_y"} and \code{"free"}.
#' @param shape.by Variable name (from \code{colData(fcs.SCE)}) for dot shaping. Default = \code{NULL}.
#' @param y.limits Numeric vector with limits for y-axis (minimum, maximum). Default = \code{NULL}.
#' @param show.stats Significances should be added to boxplots? Default = \code{TRUE}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{TRUE}.
#' @param plot.only.sig Vector indicating if only significant cell clusters should be displayed (logical element) and the P-value cutoff for selecting those ones (numerical element). Default = \code{c(FALSE, 0.05)}.
#' @param colors Vector with colors for plotting. Default = \code{NULL} (i.e., it will choose automatically a vector of colors according to \code{\link[FlowCT.v2:div.colors]{FlowCT.v2::div.colors()}}).
#' @keywords differential boxplot
#' @keywords cell clusters distributions
#' @export boxplot.cell.clustering
#' @import ggplot2
#' @importFrom data.table melt as.data.table
#' @importFrom matrixTests col_kruskalwallis
#' @importFrom stats pairwise.wilcox.test
#' @importFrom ggpubr stat_compare_means
#' @examples
#' \dontrun{
#' # option 1: show all cell clusters and return statistics
#' bx_sig <- boxplot.cell.clustering(fcs.SCE = fcs, cell.clusters = fcs$SOM_named, facet = T, 
#'     return.stats = T, 
#'     facet.free.scale = "free")
#' 
#' # option 2: show only those significant cell clusters
#' boxplot.cell.clustering(fcs.SCE = fcs, cell.clusters = fcs$SOM_named, return.stats = F, 
#'     plot.only.sig = c(T, 0.1))
#' }

boxplot.cell.clustering <- function(fcs.SCE, assay.i = "normalized", cell.clusters, condition.column = "condition", 
                                    pvalue.cutoffs = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("xxxx", "***", "**", "*", "ns")), 
                                    color.by = condition.column, geom.point = T, facet = F, facet.free.scale = "free_x", shape.by = NULL, y.limits = NULL,
                                    show.stats = T, return.stats = T, plot.only.sig = c(F, 0.05), colors = NULL){
  ## prepare tables
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SCE = fcs.SCE, assay.i,
                                                         cell.clusters = cell.clusters, count.by = "filename", plot = F)))
  
  prop_table_md <- merge(fcs.SCE@metadata$reduced_metada, prop_table, by.x = "filename", by.y = "row.names")
  
  prop_table_mdm <- melt(as.data.table(prop_table_md), id.vars = colnames(fcs.SCE@metadata$reduced_metada), 
                         variable.name = "cell_cluster", value.name = "proportion")
  
  ## statistics table
  resultskw <- col_kruskalwallis(prop_table_md[,as.character(unique(cell.clusters))], prop_table_md[,condition.column])
  KWsig <- rownames(resultskw[resultskw$pvalue < plot.only.sig[2],])
  
  if(length(unique(prop_table_md[,condition.column])) > 2){
    kw_posthoc <- list()
    for(i in KWsig) kw_posthoc[[i]] <- pairwise.wilcox.test(prop_table_md[,i], prop_table_md[,condition.column])
  }
  
  ## plotting
  if(is.null(colors)) colors <- div.colors(length(unique(prop_table_md[,condition.column])))
  
  if(plot.only.sig[1]){
    prop_table_mdm <- prop_table_mdm[prop_table_mdm$cell_cluster %in% KWsig,]
  }
  
  g <- ggplot(prop_table_mdm, aes_string("cell_cluster", "proportion", color = color.by, fill = color.by)) + 
    geom_boxplot(alpha = 0.6) + labs(x = "cell clusters", y = "Proportion") +
    theme_bw()
  
  if(!is.null(y.limits)) g <- g + scale_y_continuous(limits = c(y.limits))
  
  if(geom.point) g <- g + geom_point(aes_string("cell_cluster", "proportion", color = color.by), 
                                     alpha = 0.8, position = position_jitterdodge())
  
  if(!is.null(shape.by)) g <- g + geom_point(aes_string("cell_cluster", "proportion", color = color.by,
                                                        shape = shape.by), alpha = 0.8, position = position_jitterdodge())
  
  if(facet){
    g <- g + facet_wrap(~ cell_cluster, scales = facet.free.scale) +
                  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  }else{
    g <- g + coord_flip()
  } 
  
  if(show.stats) g <- g + stat_compare_means(label = "p.signif", symnum.args = pvalue.cutoffs)
  
  print(g + scale_fill_manual(values = colors) + scale_color_manual(values = colors))
  
  if(return.stats & length(unique(prop_table_md[,condition.column])) > 2){
    return(list(kruskal_results = resultskw, kw_posthoc_results = kw_posthoc))
  }else if(return.stats){
    return(resultskw)
  }
}
