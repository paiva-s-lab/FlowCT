#' Boxplots for identified clusters
#'
#' It draws a boxplot with cell clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}}.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters Name of column containing clusters identified through \code{\link[FlowCT:clustering.flow]{FlowCT::clustering.flow()}}.
#' @param condition Column name from the \code{colData(fcs.SCE)} object which contains condition information.
#' @param return.mode String for specifying if final resuls should be proportions ("percentage") or raw counts ("counts"). Default = \code{"percentage"}.
#' @param log.trans Logarithmic transformation of counts/percentage values?. Default = \code{FALSE}.
#' @param color.by Variable name (from \code{colData(fcs.SCE)}) for lines coloring. Default, same as \code{condition}.
#' @param facet Logical indicating if splitting boxplots by cell clusters. Default = \code{FALSE}.
#' @param facet.free.scale If \code{facet = TRUE}, string indicating how scales would be shared across all facets. Possible values: \code{"free_x"} (default), \code{"free_y"} and \code{"free"}.
#' @param y.limits Numeric vector with limits for y-axis (minimum, maximum). Default = \code{NULL}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{FALSE}.
#' @param plot.only.sig Vector indicating if only significant cell clusters should be displayed (logical element) and the P-value cutoff for selecting those ones (numerical element). Default = \code{c(FALSE, 0.05)}.
#' @param colors Vector with colors for plotting option \code{condition}.
#' @param colors Vector with colors for plotting option \code{color.by}.
#' @keywords differential boxplot
#' @keywords cell clusters distributions
#' @export boxplot.cell.clustering
#' @import ggplot2
#' @examples
#' \dontrun{
#' # option 1: show all cell clusters and return statistics
#' bx_sig <- boxplot.cell.clustering(fcs.SCE = fcs, cell.clusters = "SOM_named", facet = T,
#'     facet.free.scale = "free", return.stats = T)
#'
#' # option 2: show only those significant cell clusters
#' boxplot.cell.clustering(fcs.SCE = fcs, cell.clusters = "SOM_named",
#'     plot.only.sig = c(T, 0.1))
#' }

boxplot.cell.clustering <- function(fcs.SCE, assay.i = "normalized", cell.clusters, condition = "condition",
                                    return.mode = "percentage", log.trans = F,
                                    color.by = condition, facet = F, facet.free.scale = "free_x", y.limits,
                                    return.stats = F, plot.only.sig = c(F, 0.05), cond.colors, color.by.colors){
  ## prepare tables
  prop_table2 <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SCE, cell.clusters = cell.clusters, count.by = "filename", plot = F, assay.i = assay.i, return.mode = return.mode)))
  prop_table_md2 <- merge(distinct(as.data.frame(colData(fcs.SCE)), .data$filename, .keep_all = T) %>% select(-cell.clusters), 
                          reshape2::melt(as.matrix(prop_table2), measure.vars = as.vector(unique(fcs.SCE[[cell.clusters]])), varnames = c("filename", cell.clusters)), 
                          by = "filename")

  ## statistics
  sigs <- prop_table_md2 %>% group_by(.dots = cell.clusters) %>% rstatix::kruskal_test(value ~ condition) #.dots, to use string format
  sigs_phoc <- prop_table_md2 %>% group_by(.dots = cell.clusters) %>% rstatix::dunn_test(value ~ condition) %>% select(-c(.y.))

  ## plotting
  prop_table_md2[,cell.clusters] <- factor(prop_table_md2[,cell.clusters])
  if(missing(cond.colors)) cond.colors <- div.colors(length(unique(prop_table_md2[,condition])))
  if(missing(color.by.colors)) color.by.colors <- div.colors(length(unique(prop_table_md2[,color.by])))

  if(plot.only.sig[1]){
    prop_table_md2 <- prop_table_md2[prop_table_md2[,cell.clusters] %in% sigs[cell.clusters][sigs[,"p"] <= plot.only.sig[2]],]
  }

  g <- ggplot(prop_table_md2, aes_string(cell.clusters, "value", fill = condition)) +
    geom_boxplot(aes_string(color = condition), alpha = 0.6) + 
    # geom_point(aes_string(cell.clusters, "value", color = color.by), position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.6)) +
    # geom_jitter(aes_string(cell.clusters, "value", color = color.by), width = 0.2, height = 0.6) +
    scale_fill_manual(values = cond.colors, guide = F) + scale_color_manual(values = c(cond.colors, color.by.colors)) +
    labs(x = "cell clusters", y = "Proportion") + theme_bw() + theme(legend.position = "bottom") +
    ggpubr::stat_compare_means(label = "p.signif")

  if(return.mode == "percentage") g <- g + geom_jitter(aes_string(cell.clusters, "value", color = color.by), width = 0.2, height = 0.6*10) else
    g <- g + geom_jitter(aes_string(cell.clusters, "value", color = color.by), width = 0.2, height = 0.6)

  if(!missing(y.limits)) g <- g + scale_y_continuous(limits = c(y.limits))

  if(facet){
    g <- g + facet_wrap(as.formula(paste("~", cell.clusters)), scales = facet.free.scale) +
                  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  }else{
    g <- g + coord_flip()
  }

  if(log.trans){
    g <- g + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                        labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
             labs(x = cell.clusters, y = paste(return.mode, "of cells (log10-transf.)\n"), color = "")
  }else{
     g <- g + labs(x = cell.clusters, y = paste(return.mode, "of cells\n"), color = "") 
  }

  if(return.stats & length(unique(prop_table_md2[,condition])) > 2){
    print(g)
    return(list(kruskal = sigs, kw_posthoc = sigs_phoc, plot = g))
  }else if(return.stats){
    print(g)
    return(list(kruskal_results = sigs, plot = g))
  }else{
    return(g)
  }
}
