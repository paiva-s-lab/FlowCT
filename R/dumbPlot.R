#' Dumbbell plot with clusters
#'
#' It draws a Dumbbell plot according condition for each cell cluster identified.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param condition.column Column name from the \code{colData(fcs.SCE)} object which contains condition information.
#' @param psig.cutoff P-value cutoff. Default = \code{0.05}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{FALSE}.
#' @param colors Vector with colors for plotting (same length as different conditions within the experiment).
#' @keywords differential dotplot
#' @keywords Dumbbell plot
#' @keywords longitudinal dotplot
#' @export
#' @import ggplot2
#' @import dplyr
#' @examples
#' \dontrun{
#' diffDots.cell.clustering(fcs.SCE = fcs, cell.clusters = "SOM_named", return.stats = T)
#' }

dumbPlot <- function(fcs.SCE, assay.i = "normalized", cell.clusters, condition, psig.cutoff = 0.05, return.stats = F, colors){
  if(length(unique(prop_table_ms[,condition]))) stop("Not posible compare only a condition.", call. = F)
  
  ## prepare tables
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SCE, cell.clusters = cell.clusters, count.by = condition, plot = F, assay.i = assay.i, return.mode = "percentage")))
  prop_table_m <- reshape2::melt(as.matrix(prop_table), measure.vars = as.vector(unique(fcs.SCE[[cell.clusters]])), varnames = c(condition, cell.clusters))

  prop_table2 <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SCE, cell.clusters = cell.clusters, count.by = "filename", plot = F, assay.i = assay.i, return.mode = "counts")))
  prop_table_md2 <- merge(distinct(as.data.frame(colData(fcs.SCE)), .data$filename, .keep_all = T) %>% select(-cell.clusters), 
                          reshape2::melt(as.matrix(prop_table2), measure.vars = as.vector(unique(fcs1[[cell.clusters]])), varnames = c("filename", cell.clusters)), 
                          by = "filename")

  ## statistics
  sigs <- prop_table_md2 %>% group_by(.dots = cell.clusters) %>% rstatix::kruskal_test(value ~ condition) #.dots, to use string format
  prop_table_ms <- merge(prop_table_m, sigs[,c(cell.clusters, "p")], by = cell.clusters)
  prop_table_ms$sigcol <- ifelse(prop_table_ms$p <= psig.cutoff, "sig", "ns")

  sigs_phoc <- prop_table_md2 %>% group_by(.dots = cell.clusters) %>% rstatix::dunn_test(value ~ condition) %>% select(-c(.y.))

  ## plotting
  if(missing(colors)) colors <- div.colors(length(unique(prop_table_ms[,condition])))

  g <- ggplot(prop_table_ms, aes_string(x = "value", y = "SOM.k10", fill = "condition", color = "sigcol")) +
          geom_line(aes_string(group = "SOM.k10"), size = 1) +
          geom_point(size = 4, shape = 21, color = "white") +
          scale_color_manual(values = c("gray63", "brown1"), labels = c("no sig.", "sig.")) +
          scale_fill_manual(values = colors, labels = unique(prop_table_ms[,"condition"])) +
          scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
          theme(panel.background = element_blank(), legend.key = element_blank(),
                panel.grid.major.x = element_line(colour = "gray73", linetype = "dashed", size = 0.3),
                panel.border = element_rect(color = "gray73", fill = NA),
                axis.text = element_text(color = "gray20", face = "bold"), axis.title = element_text(face = "bold")) +
          xlab("\n Population % (log10-transformed)") + ylab("Cell clusters\n") + labs(color = "")

  if(return.stats & length(colors) > 2){
    return(list(kruskal = sigs, kw_posthoc = sigs_phoc, plot = g))
  }else if(return.stats){
    return(list(kruskal_results = sigs, plot = g))
  }else{
    return(g)
  }
}
