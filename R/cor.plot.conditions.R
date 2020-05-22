#' cor.plot.conditions
#'
#' It draws a correlation plot (between two specified conditions within a \code{fcs.SCE} object.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT]{fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT]{fsom.clustering()}} (and, normaly, later renamed).
#' @param condition.column Column name from the \code{colData(fcs.SCE)} object which contains condition information.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{FALSE}.
#' @keywords correlation plot
#' @keywords corr
#' @export
#' @importFrom Hmisc rcorr
#' @importFrom corrplot corrplot
#' @examples
#' \dontrun{
#' corplot.conditions(fcs.SCE = fcs, cell.clusters = fcs$SOM_named, condition.column = "condition")
#' }

corplot.conditions <- function(fcs.SCE, assay.i = "normalized", cell.clusters, condition.column, return.stats = F){
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SCE = fcs.SCE, cell.clusters = cell.clusters, count.by = "filename", plot = F)))
  
  prop_table_md <- merge(fcs.SCE@metadata$reduced_metada, prop_table, by.x = "filename", by.y = "row.names")
  
  conditions <- as.character(unique(prop_table_md[,condition.column]))
  if(length(conditions) != 2) stop("corplot is onlyprepared for comparing ONLY two conditions")
  
  dataset1 <- prop_table_md[prop_table_md[,condition.column] == conditions[1],colnames(prop_table)]
  colnames(dataset1) <- paste0(colnames(dataset1),  ":", conditions[1])
  rownames(dataset1) <- NULL
  dataset2 <- prop_table_md[prop_table_md[,condition.column] == conditions[2],colnames(prop_table)]
  colnames(dataset2) <- paste0(colnames(dataset2), ":", conditions[2])
  rownames(dataset2) <- NULL
  
  dataset <- merge(dataset1, dataset2, by = "row.names")[-1]
  corr <- rcorr(as.matrix(dataset))
  
  # face conditions
  corr_r <- as.matrix(corr$r[grepl(conditions[1], rownames(corr$r)), grepl(conditions[2], colnames(corr$r))]) 
  pval <- as.matrix(corr$P[grepl(conditions[1], rownames(corr$P)), grepl(conditions[2], colnames(corr$P))])
  
  corrplot(corr_r, order = "hclust", p.mat = pval, insig = "label_sig",
                           sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white",
                           tl.col="black", tl.cex=.7, tl.offset=0.5, tl.srt=45, 
                           col = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                        "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                        "#4393C3", "#2166AC", "#053061")))(100))
  if(return.stats) return(corr_r)
}

