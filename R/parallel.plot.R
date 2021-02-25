#' Longitudinal differential dotplot
#'
#' It draws a differential dot plot (longitudinaly) according condition for each cell cluster identified.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param condition Column name from the \code{colData(fcs.SCE)} object which contains condition information.
#' @param psig.cutoff P-value cutoff. Default = \code{0.05}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{FALSE}.
#' @param colors Vector with colors for plotting.
#' @param return.mode String for specifying if final resuls should be proportions ("percentage") or raw counts ("counts"). Default = \code{"percentage"}.
#' @param hide.nosig Logical indicating whether hiding non-significal cell populations. Default = \code{FALSE}.
#' @param log.trans Logarithmic transformation of counts/percentage values?. Default = \code{FALSE}.
#' @param size Point size. Default = \code{3}.
#' @param labels.pos Position for cell clusters labelling, user should indicate the numeric position (i.e., 1 for first condition, 2 for second, and so on). By default, last condition will be used.
#' @keywords differential parallel dotplot
#' @export
#' @import ggplot2
#' @import dplyr
#' @examples
#' \dontrun{
#' parallel.plot(fcs.SCE = fcs, cell.clusters = "SOM_named", condition = "condition")
#' }

parallel.plot <- function(fcs.SCE, assay.i = "normalized", cell.clusters, condition, psig.cutoff = 0.05, return.stats = F, colors, return.mode = "percentage", hide.nosig = F, log.trans = F, size = 3, labels.pos){
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package \"ggrepel\" needed for this function to work. Please install it.", call. = FALSE)

  conditions <- unique(fcs.SCE[[condition]])
  if(length(conditions) == 1) stop("Not posible compare only a condition.", call. = F)
  
  ## prepare tables
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SCE, cell.clusters = cell.clusters, count.by = condition, plot = F, assay.i = assay.i, return.mode = return.mode)))
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
  if(missing(colors)) colors <- div.colors(length(unique(fcs.SCE[[cell.clusters]])))
  prop_table_ms[,cell.clusters] <- factor(prop_table_ms[,cell.clusters])

  g <- ggplot(prop_table_ms, aes_string(x = condition, y = "value", fill = cell.clusters)) +
          scale_fill_manual(values = colors, na.value  = "gray63") +
          theme(panel.background = element_blank(), axis.line = element_line(color = "black"))

    if(hide.nosig){
      g <- g + geom_line(aes_string(group = cell.clusters, color = "sigcol"), size = 1) +
           geom_point(data = subset(prop_table_ms, prop_table_ms$sigcol == "ns"), size = size, shape = 21, color = "white", fill = "gray63") + 
           geom_point(data = subset(prop_table_ms, prop_table_ms$sigcol == "sig"), size = size, shape = 21, color = "white") 
    }else{
      g <- g + geom_line(aes_string(group = cell.clusters, color = "sigcol")) +
            geom_point(aes_string(fill = cell.clusters), size = size, shape = 21, color = "white")
    }

  if(sum(prop_table_ms$sigcol == "sig") == 0){
    g <- g + scale_colour_manual(values = c("gray63"), labels = c("no sig."))  
  }else{
    g <- g + scale_colour_manual(values = c("gray63", "brown1"), labels = c("no sig.", "sig."))
  }

  if(log.trans){
    g <- g + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                        labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
             labs(x = "\nCondition", y = paste(return.mode, "of cells (log10-transf.)\n"), color = "", fill = "Cell clusters")
  }else{
     g <- g + labs(x = "\nCondition", y = paste(return.mode, "of cells\n"), color = "", fill = "Cell clusters") 
  }

  if(missing(labels.pos)) labels.pos <- length(conditions)
  g <- g + ggrepel::geom_label_repel(data = subset(prop_table_ms, prop_table_ms$condition == conditions[labels.pos]), 
                            aes_string(x = labels.pos, y = "value", label = cell.clusters, fill = cell.clusters),
                            nudge_x = 0.2, show.legend = F)


  if(return.stats & length(conditions) > 2){
    print(g)
    return(list(kruskal = sigs, kw_posthoc = sigs_phoc, plot = g))
  }else if(return.stats){
    print(g)
    return(list(kruskal_results = sigs, plot = g))
  }else{
    return(g)
  }
}
