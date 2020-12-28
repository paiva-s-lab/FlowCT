#' Longitudinal differential dotplot
#'
#' It draws a differential dot plot (longitudinaly) according condition for each cell cluster identified.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param condition.column Column name from the \code{colData(fcs.SCE)} object which contains condition information. De
#' @param psig.cutoff P-value cutoff. Default = \code{0.05}.
#' @param colors Vector with colors for plotting. Default = \code{NULL} (i.e., it will choose automatically a vector of colors according to \code{\link[FlowCT:div.colors]{FlowCT::div.colors()}}).
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{TRUE}.
#' @param hide.nosig Logical indicating whether hiding non-significal cell populations. Default = \code{FALSE}.
#' @keywords differential dotplot
#' @keywords Dumbbell plot
#' @keywords longitudinal dotplot
#' @export
#' @import ggplot2
#' @importFrom data.table melt as.data.table
#' @importFrom stats aggregate median pairwise.wilcox.test ave
#' @importFrom matrixTests col_kruskalwallis
#' @examples
#' \dontrun{
#' diffDots.cell.clustering(fcs.SCE = fcs, cell.clusters = fcs$SOM_named, return.stats = F)
#' }

diffdots.cell.clustering <- function(fcs.SCE, assay.i = "normalized", cell.clusters, condition.column, psig.cutoff = 0.05, return.stats = T, colors = NULL, return.mode = "percentage", hide.nosig = F){
  ## prepare tables
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SCE, cell.clusters, count.by = "filename", plot = F, assay.i = assay.i, return.mode = return.mode)))

  prop_table_md <- merge(fcs.SCE@metadata$reduced_metada, prop_table, by.x = "filename", by.y = "row.names")

  dfm <- as.data.frame(melt(as.data.table(prop_table_md), measure.vars = as.vector(unique(fcs.SCE[[cell.clusters]]))))
  dfma <- aggregate(dfm$value ~ dfm$variable + dfm[,condition.column], FUN = median)
  colnames(dfma) <- c("variable", "condition", "value")
  dfma <- transform(dfma, pct = log(ave(dfma$value, dfma$condition, FUN = function(x) x/sum(x)*100))) #transform to percentaje

  ## keep original factor order for plotting
  if(class(fcs.SCE[[condition.column]]) == "factor") dfma$condition <- factor(dfma$condition, levels = levels(fcs.SCE[[condition.column]]))
  if(class(fcs.SCE[[cell.clusters]]) == "factor") dfma$variable <- factor(dfma$variable, levels = levels(fcs.SCE[[cell.clusters]]))

  conditions <- unique(dfma[,condition.column])

  ## statistics table
  resultskw <- col_kruskalwallis(prop_table_md[,as.character(unique(fcs.SCE[[cell.clusters]]))], prop_table_md[,condition.column])
  KWsig <- rownames(resultskw[resultskw$pvalue < psig.cutoff,])
  dfma$sig <- ifelse(dfma$variable %in% KWsig, "1", "0")

  if(length(unique(prop_table_md[,condition.column])) > 2){
    kw_posthoc <- list()
    for(i in KWsig) kw_posthoc[[i]] <- pairwise.wilcox.test(prop_table_md[,i], prop_table_md[,condition.column])
  }

  ## plotting
  if(is.null(colors)) colors <- div.colors(length(unique(fcs.SCE[[cell.clusters]])))
  g <- ggplot(dfma, aes_string(x = "condition", y = "pct", fill = "variable")) +
          scale_colour_manual(values = c("gray63", "brown1"), labels = c("no sig.", "sig.")) +
          scale_fill_manual(values = colors, na.value  = "gray63") +
          theme(panel.background = element_blank(), axis.line = element_line(color = "black")) +
          labs(x = "\nCondition", y = "% of cells (log-transf.)\n", color = "", fill = "Cell clusters")

  if(!hide.nosig){
     g <- g + geom_line(aes_string(group = "variable", color = "sig"), size = 1) +
                          geom_point(data = subset(dfma, dfma$sig == 0), size = 3, shape = 21, color = "gray63", fill = "gray63") + 
                          geom_point(data = subset(dfma, dfma$sig == 1), size = 3, shape = 21, color = "gray63")
  }else{
     g <- g + geom_line(data = subset(dfma, dfma$sig == 1), aes_string(group = "variable", color = "sig"), size = 1) +
              geom_point(data = subset(dfma, dfma$sig == 1), size = 3, shape = 21, color = "gray63")
  }

  if(sum(dfma$sig == 1) != 0){
    if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package \"ggrepel\" needed for this function to work. Please install it.", call. = FALSE)
    g <- g + ggrepel::geom_label_repel(data = subset(dfma, dfma$sig == 1 & dfma$condition == conditions[length(conditions)]),
                            aes_string(x = 1, y = "pct", label = "variable", fill = "variable"),
                            nudge_x = -0.1, show.legend = F)
  }


  if(return.stats & length(unique(prop_table_md[,condition.column])) > 2){
    return(list(kruskal_results = resultskw, kw_posthoc_results = kw_posthoc, plot = g))
  }else if(return.stats){
    return(list(kruskal_results = resultskw, plot = g))
  }else{
    return(g)
  }
}
