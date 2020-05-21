# 'dumbPlot.cell.clustering
#'
#' It draws a Dumbbell plot according condition for each cell cluster identified.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param condition.column Column name from the \code{colData(fcs.SCE)} object which contains condition information. De
#' @param psig.cutoff P-value cutoff. Default = \code{0.05}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{TRUE}.
#' @keywords differential dotplot
#' @keywords Dumbbell plot
#' @keywords longitudinal dotplot
#' @export 
#' @import ggplot2
#' @importFrom stats aggregate median pairwise.wilcox.test ave
#' @importFrom data.table melt
#' @importFrom matrixTests col_kruskalwallis
#' @examples
#' \dontrun{
#' diffDots.cell.clustering(fcs.SCE = fcs.SCE, cell.clusters = fcs.SCE$SOM_named, return.stats = F)
#' }

dumbPlot.cell.clustering <- function(fcs.SCE, assay.i = "normalized", cell.clusters, condition.column, psig.cutoff = 0.05, return.stats = T){
  ## prepare tables
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SCE, cell.clusters, count.by = "filename", plot = F, assay.i = assay.i)))
  
  prop_table_md <- merge(fcs.SCE@metadata$reduced_metada, prop_table, by.x = "filename", by.y = "row.names")
  
  dfm <- melt(prop_table_md, measure.vars = unique(cell.clusters))
  dfma <- aggregate(dfm$value ~ dfm$variable + dfm[,condition.column], FUN = median)
  colnames(dfma) <- c("variable", "condition", "value")
  dfma <- transform(dfma, pct = log(ave(dfma$value, dfma[,condition.column], FUN = function(x) x/sum(x)*100))) #transform to percentaje
  
  ## statistics table
  resultskw <- col_kruskalwallis(prop_table_md[,as.character(unique(cell.clusters))], prop_table_md[,condition.column])
  KWsig <- rownames(resultskw[resultskw$pvalue < psig.cutoff,])
  dfma$sig <- ifelse(dfma$variable %in% KWsig, "1", "0")
  
  if(length(unique(prop_table_md[,condition.column])) > 2){
    kw_posthoc <- list()
    for(i in KWsig) kw_posthoc[[i]] <- pairwise.wilcox.test(prop_table_md[,i], prop_table_md[,condition.column])
  }
  
  # beta: calculating post-hoc for multiple conditions and draw colored lines between each condition and not general
  # KW_ph_sig <- c()
  # if(length(unique(prop_table_md[,"condition"])) > 2){
  #   for(i in names(kw_posthoc)){
  #     aux_stats <- melt(kw_posthoc[[i]]$p.value)
  #     aux_stats <- aux_stats[aux_stats$value < 0.1,]
  #   }
  # }
  
  ## plotting
  print(ggplot(dfma, aes_string(x = "pct", y = "variable", fill = condition.column, color = "sig")) + 
          geom_line(aes_string(group = "variable"), size = 1) +
          scale_color_manual(values = c("gray63", "brown1"), labels = c("no sig.", "sig.")) +
          geom_point(size = 4, shape = 21, color = "white") +
          scale_fill_manual(values = div.colors(length(unique(prop_table_md[,condition.column])), set.seed = 3), 
                            labels = unique(prop_table_md[,condition.column])) +
          # guides(color = guide_legend(ncol = 1)) + #display legend in one-column format
          theme(panel.background = element_blank(), legend.key = element_blank(),
                panel.grid.major.x = element_line(colour = "gray73", linetype = "dashed", size = 0.3),
                panel.border = element_rect(color = "gray73", fill = NA), 
                axis.text = element_text(color = "gray20", face = "bold"), axis.title = element_text(face = "bold")) +
          xlab("\n Population % (log-transformed)") + ylab("Cell clusters\n") + labs(color = ""))
  
  if(return.stats & length(unique(prop_table_md[,condition.column])) > 2){
    return(list(kruskal_results = resultskw, kw_posthoc_results = kw_posthoc))
  }else if(return.stats){
    return(resultskw)
  }
}
